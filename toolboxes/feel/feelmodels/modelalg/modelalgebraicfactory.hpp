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
        typedef std::shared_ptr<model_type> model_ptrtype;
        typedef std::weak_ptr<model_type> model_weakptrtype;

        using index_type = typename model_type::index_type;
        typedef typename model_type::value_type value_type;
        typedef typename model_type::backend_type backend_type;
        typedef typename model_type::backend_ptrtype backend_ptrtype;

        typedef GraphCSR graph_type;
        typedef std::shared_ptr<graph_type> graph_ptrtype;

        typedef typename model_type::vector_ptrtype vector_ptrtype;
        typedef typename model_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
        typedef typename sparse_matrix_ptrtype::element_type sparse_matrix_type;

        typedef Preconditioner<value_type> preconditioner_type;
        typedef std::shared_ptr<preconditioner_type> preconditioner_ptrtype;

        typedef typename backend_type::solve_return_type solve_return_type;
        typedef typename backend_type::nl_solve_return_type nl_solve_return_type;

        typedef model_type::indexsplit_type indexsplit_type;
        typedef model_type::indexsplit_ptrtype indexsplit_ptrtype;


        typedef boost::function<void ( ModelAlgebraic::DataUpdateLinear& )> function_assembly_linear_type;
        typedef boost::function<void ( ModelAlgebraic::DataUpdateJacobian& )> function_assembly_jacobian_type;
        typedef boost::function<void ( ModelAlgebraic::DataUpdateResidual& )> function_assembly_residual_type;
        typedef boost::function<void ( ModelAlgebraic::DataNewtonInitialGuess& )> function_newton_initial_guess_type;

        typedef typename backend_type::pre_solve_type pre_solve_type;
        typedef typename backend_type::post_solve_type post_solve_type;
        typedef typename backend_type::update_nlsolve_type update_nlsolve_type;

        //---------------------------------------------------------------------------------------------------------------//
        //---------------------------------------------------------------------------------------------------------------//
        //---------------------------------------------------------------------------------------------------------------//

        ModelAlgebraicFactory( std::string const& prefix, po::variables_map const& vm = Environment::vm() );
        ModelAlgebraicFactory( model_ptrtype const& model, backend_ptrtype const& backend );

        ModelAlgebraicFactory(model_ptrtype const& model,
                              backend_ptrtype const& backend,
                              graph_ptrtype const& graph,
                              indexsplit_ptrtype const& indexSplit );

        //---------------------------------------------------------------------------------------------------------------//
        //---------------------------------------------------------------------------------------------------------------//

        void init( model_ptrtype const& model, backend_ptrtype const& backend );
        void init( model_ptrtype const& model, backend_ptrtype const& backend,
                   graph_ptrtype const& graph, indexsplit_ptrtype const& indexSplit );
        void init( backend_ptrtype const& backend, graph_ptrtype const& graph, indexsplit_ptrtype const& indexSplit );

        void initExplictPartOfSolution();
        vector_ptrtype explictPartOfSolution() { return M_explictPartOfSolution; }

        bool useSolverPtAP() const { return M_useSolverPtAP; }
        void initSolverPtAP( sparse_matrix_ptrtype matP, sparse_matrix_ptrtype matQ );
        bool hasInitSolverPtAP() const { return M_solverPtAP_backend? true : false; }
        void solverPtAP_setDofEliminationIds( std::set<index_type> const& dofId ) { M_solverPtAP_dofEliminationIds = dofId; }
        sparse_matrix_ptrtype solverPtAP_matrixP() const { return M_solverPtAP_matP; }

        //---------------------------------------------------------------------------------------------------------------//

        void
        rebuildMatrixVector( graph_ptrtype const& graph,
                             indexsplit_ptrtype const& indexSplit);

        void reset(backend_ptrtype __backend,
                   graph_ptrtype const& graph,
                   indexsplit_ptrtype const& indexSplit);

        std::shared_ptr<std::ostringstream> getInfo() const;
        void printInfo() const;

        void updateInformationObject( nl::json & p ) const;
        static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

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

        //! return data information available in assembly process
        ModelAlgebraic::DataUpdateBase const& dataInfos() const { return M_dataInfos; }

        //! return data information available in assembly process
        ModelAlgebraic::DataUpdateBase & dataInfos() { return M_dataInfos; }

        void
        attachNullSpace( NullSpace<value_type> const& nullSpace );
        void
        attachNearNullSpace( NullSpace<value_type> const& nearNullSpace );
        void
        attachNearNullSpace( int k, NullSpace<value_type> const& nearNullSpace, std::string const& prefix );
        void
        attachNearNullSpace( int k, NullSpace<value_type> const& nearNullSpace );

        //! attach a sparse matrix to the precondtioner
        void attachAuxiliarySparseMatrix( std::string const& key,sparse_matrix_ptrtype const& mat );
        //! return true if a sparse matrix has been attached
        bool hasAuxiliarySparseMatrix( std::string const& key ) const;
        //! return  a sparse matrix attached
        sparse_matrix_ptrtype const& auxiliarySparseMatrix( std::string const& key ) const;

        //! attach operator PCD to the precondtioner
        void attachOperatorPCD( std::string const& key, typename preconditioner_type::operator_pcdbase_ptrtype const& opPCD );

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

        void preSolveNewton( vector_ptrtype rhs, vector_ptrtype sol ) const;
        void postSolveNewton( vector_ptrtype rhs, vector_ptrtype sol ) const;
        void preSolvePicard( vector_ptrtype rhs, vector_ptrtype sol ) const;
        void postSolvePicard( vector_ptrtype rhs, vector_ptrtype sol ) const;
        void preSolveLinear( vector_ptrtype rhs, vector_ptrtype sol ) const;
        void postSolveLinear( vector_ptrtype rhs, vector_ptrtype sol ) const;


        void rebuildCstJacobian( vector_ptrtype U );
        void rebuildCstLinearPDE( vector_ptrtype U );


        //! apply assembly of linear operators rhs and lhs (can be usefull for an external use)
        void applyAssemblyLinear(const vector_ptrtype& U, sparse_matrix_ptrtype& lhs, vector_ptrtype& rhs,
                                 std::vector<std::string> const& infos = std::vector<std::string>(),
                                 bool applyDofElimination = true ) const;
        void applyAssemblyLinear( ModelAlgebraic::DataUpdateLinear & dataLinear, bool applyDofElimination = true ) const;

        void evaluateResidual( const vector_ptrtype& U, vector_ptrtype& R,
                               std::vector<std::string> const& infos = std::vector<std::string>(),
                               bool applyDofElimination = true ) const;
        void evaluateResidual( ModelAlgebraic::DataUpdateResidual & dataResidual, bool applyDofElimination = true ) const;
        //---------------------------------------------------------------------------------------------------------------//
        //---------------------------------------------------------------------------------------------------------------//
        //---------------------------------------------------------------------------------------------------------------//

        void setFunctionLinearAssembly( function_assembly_linear_type const& func ) { M_functionLinearAssembly = func; }
        void setFunctionLinearDofElimination( function_assembly_linear_type const& func ) { M_functionLinearDofElimination = func; }
        void setFunctionNewtonInitialGuess( function_newton_initial_guess_type const& func ) { M_functionNewtonInitialGuess = func; }
        void setFunctionJacobianAssembly( function_assembly_jacobian_type const& func ) { M_functionJacobianAssembly = func; }
        void setFunctionResidualAssembly( function_assembly_residual_type const& func ) { M_functionResidualAssembly = func; }

        void addFunctionLinearAssembly( function_assembly_linear_type const& func, std::string const& key = "" );
        void addFunctionLinearDofElimination( function_assembly_linear_type const& func, std::string const& key = "" );
        void addFunctionLinearPostAssembly( function_assembly_linear_type const& func, std::string const& key = "" );
        void addFunctionNewtonInitialGuess( function_newton_initial_guess_type const& func, std::string const& key = "" );
        void addFunctionJacobianAssembly( function_assembly_jacobian_type const& func, std::string const& key = "" );
        void addFunctionResidualAssembly( function_assembly_residual_type const& func, std::string const& key = "" );
        void addFunctionJacobianDofElimination( function_assembly_jacobian_type const& func, std::string const& key = "" );
        void addFunctionResidualDofElimination( function_assembly_residual_type const& func, std::string const& key = "" );
        void addFunctionJacobianPostAssembly( function_assembly_jacobian_type const& func, std::string const& key = "" );
        void addFunctionResidualPostAssembly( function_assembly_residual_type const& func, std::string const& key = "" );

        void addVectorLinearRhsAssembly( vector_ptrtype const& vec, double scaling = 1.0, std::string const& key = "", bool cstPart = false );
        void addVectorResidualAssembly( vector_ptrtype const& vec, double scaling = 1.0, std::string const& key = "", bool cstPart = false );
        void setActivationAddVectorLinearRhsAssembly( std::string const& key, bool b );
        void setActivationAddVectorResidualAssembly( std::string const& key, bool b );

        void updateNewtonIteration( int step, vector_ptrtype residual, vector_ptrtype sol, typename backend_type::solvernonlinear_type::UpdateIterationData const& data );
        void updatePicardIteration( int step, vector_ptrtype sol );
    private :

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

        vector_ptrtype M_explictPartOfSolution;
        vector_ptrtype M_contributionsExplictPartOfSolutionWithNewton;

        bool M_useSolverPtAP;
        sparse_matrix_ptrtype M_solverPtAP_matP;
        sparse_matrix_ptrtype M_solverPtAP_matQ; // operator from natural basis to PtAP basis
        sparse_matrix_ptrtype M_solverPtAP_matPtAP;
        vector_ptrtype M_solverPtAP_PtF;
        vector_ptrtype M_solverPtAP_solution;
        vector_ptrtype M_solverPtAP_Psolution;
        preconditioner_ptrtype M_solverPtAP_prec;
        backend_ptrtype M_solverPtAP_backend;
        std::optional<std::set<index_type>> M_solverPtAP_dofEliminationIds;

        double M_dofElimination_valueOnDiagonal;
        Feel::Context M_dofElimination_strategy;

        bool M_hasBuildLinearJacobian;
        bool M_hasBuildResidualCst;
        bool M_hasBuildLinearSystemCst;

        ModelAlgebraic::DataUpdateBase M_dataInfos;

        function_assembly_linear_type M_functionLinearAssembly;
        function_assembly_linear_type M_functionLinearDofElimination;
        function_newton_initial_guess_type M_functionNewtonInitialGuess;
        function_assembly_jacobian_type M_functionJacobianAssembly;
        function_assembly_residual_type M_functionResidualAssembly;

        std::map<std::string,function_assembly_linear_type> M_addFunctionLinearAssembly;
        std::map<std::string,function_assembly_linear_type> M_addFunctionLinearDofElimination;
        std::map<std::string,function_assembly_linear_type> M_addFunctionLinearPostAssembly;
        std::map<std::string,function_newton_initial_guess_type> M_addFunctionNewtonInitialGuess;
        std::map<std::string,function_assembly_jacobian_type> M_addFunctionJacobianAssembly;
        std::map<std::string,function_assembly_residual_type> M_addFunctionResidualAssembly;
        std::map<std::string,function_assembly_jacobian_type> M_addFunctionJacobianDofElimination;
        std::map<std::string,function_assembly_residual_type> M_addFunctionResidualDofElimination;
        std::map<std::string,function_assembly_jacobian_type> M_addFunctionJacobianPostAssembly;
        std::map<std::string,function_assembly_residual_type> M_addFunctionResidualPostAssembly;

        // ( key -> ( vector,scaling, cstPart?, activated? ) )
        std::map<std::string, std::tuple<vector_ptrtype,double,bool,bool>> M_addVectorLinearRhsAssembly;
        std::map<std::string, std::tuple<vector_ptrtype,double,bool,bool>> M_addVectorResidualAssembly;

        bool M_usePseudoTransientContinuation;
        std::string M_pseudoTransientContinuationEvolutionMethod;
        std::vector<std::pair<double,double> > M_pseudoTransientContinuationDeltaAndResidual;
        double M_pseudoTransientContinuationDelta0, M_pseudoTransientContinuationDeltaMax;
        std::string M_pseudoTransientContinuationSerVariant;
        vector_ptrtype M_pseudoTransientContinuationPreviousSolution;
        double M_pseudoTransientContinuationExpurThresholdHigh,M_pseudoTransientContinuationExpurThresholdLow;
        double M_pseudoTransientContinuationExpurBetaHigh, M_pseudoTransientContinuationExpurBetaLow;
    };


} // namespace FeelModels
} // namespace Feel


#endif //FEELPP_MODELSALGEBRAICFACTORY_HPP
