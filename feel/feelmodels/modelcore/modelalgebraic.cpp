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
 \file modelalgebraic.cpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2012-01-19
 */

#include <feel/feelmodels/modelcore/modelalgebraic.hpp>

namespace Feel
{

namespace FeelModels
{

ModelAlgebraic::ModelAlgebraic( std::string _theprefix,
                                WorldComm const& _worldComm,
                                std::string const& subPrefix,
                                std::string const& rootRepository )
    : super_type( _theprefix, _worldComm, subPrefix, rootRepository ),
      M_verboseSolverTimer( boption( _name = "verbose_solvertimer", _prefix = this->prefix() ) ),
      M_verboseSolverTimerAllProc( boption( _name = "verbose_solvertimer_allproc", _prefix = this->prefix() ) ),
      M_rebuildCstPartInLinearSystem( boption( _name = "linearsystem-cst-update", _prefix = this->prefix() ) ),
      M_useLinearJacobianInResidual( boption( _name = "residual-uselinearjac", _prefix = this->prefix() ) ),
      M_rebuildLinearPartInJacobian( boption( _name = "jacobian-linear-update", _prefix = this->prefix() ) ),
      M_rebuildCstPartInResidual( true ), // not an option (just opitmisation with semi-implicit)
      M_useCstMatrix( boption( _name = "use-cst-matrix", _prefix = this->prefix() ) ),
      M_useCstVector( boption( _name = "use-cst-vector", _prefix = this->prefix() ) ),
      M_needToRebuildCstPart( false ),
      M_errorIfSolverNotConverged( boption( _name = "error-if-solver-not-converged", _prefix = this->prefix() ) ),
      M_printGraph( boption( _name = "graph-print-python", _prefix = this->prefix() ) )
{

    //-----------------------------------------------------------------------//
    //-----------------------------------------------------------------------//
    if ( Environment::vm().count( prefixvm( this->prefix(), "graph-print-python-filename" ) ) )
        M_printGraphFileName = Environment::vm()[prefixvm( this->prefix(), "graph-print-python-filename" )].as<std::string>();
    else
        M_printGraphFileName = this->prefix() + ".graphPython.py";
}
ModelAlgebraic::~ModelAlgebraic()
{
}

// verbose
bool ModelAlgebraic::verboseSolverTimer() const
{
    return M_verboseSolverTimer;
}
bool ModelAlgebraic::verboseSolverTimerAllProc() const
{
    return M_verboseSolverTimerAllProc;
}
// do rebuild cst part in linear/jacobian or use jac for residual
bool ModelAlgebraic::rebuildCstPartInLinearSystem() const
{
    return M_rebuildCstPartInLinearSystem;
}
void ModelAlgebraic::setRebuildCstPartInLinearSystem( bool b )
{
    M_rebuildCstPartInLinearSystem = b;
}
bool ModelAlgebraic::useLinearJacobianInResidual() const
{
    return M_useLinearJacobianInResidual;
}
void ModelAlgebraic::setUseLinearJacobianInResidual( bool b )
{
    M_useLinearJacobianInResidual = b;
}
bool ModelAlgebraic::rebuildLinearPartInJacobian() const
{
    return M_rebuildLinearPartInJacobian;
}
void ModelAlgebraic::setRebuildLinearPartInJacobian( bool b )
{
    M_rebuildLinearPartInJacobian = b;
}
// a utiliser avec precaution!!!
bool ModelAlgebraic::rebuildCstPartInResidual() const
{
    return M_rebuildCstPartInResidual;
}
void ModelAlgebraic::setRebuildCstPartInResidual( bool b )
{
    M_rebuildCstPartInResidual = b;
}
// define an other matrix/vector to store the cst part
bool ModelAlgebraic::useCstMatrix() const
{
    return M_useCstMatrix;
}
void ModelAlgebraic::setUseCstMatrix( bool b )
{
    M_useCstMatrix = b;
}
bool ModelAlgebraic::useCstVector() const
{
    return M_useCstVector;
}
void ModelAlgebraic::setUseCstVector( bool b )
{
    M_useCstVector = b;
}
// allow to rebuild cst part (once at next solve) if some parameters (model,time mode,..) change
bool ModelAlgebraic::needToRebuildCstPart() const
{
    return M_needToRebuildCstPart;
}
void ModelAlgebraic::setNeedToRebuildCstPart( bool b )
{
    M_needToRebuildCstPart = b;
}
// an option
bool ModelAlgebraic::errorIfSolverNotConverged() const
{
    return M_errorIfSolverNotConverged;
}
void ModelAlgebraic::setErrorIfSolverNotConverged( bool b )
{
    M_errorIfSolverNotConverged = b;
}
// save a python script to view graph
bool ModelAlgebraic::printGraph() const
{
    return M_printGraph;
}
void ModelAlgebraic::setPrintGraph( bool b )
{
    M_printGraph = b;
}
std::string
ModelAlgebraic::printGraphFileName() const
{
    return M_printGraphFileName;
}
void ModelAlgebraic::setPrintGraphFileName( std::string s )
{
    M_printGraphFileName = s;
}

/**
 * return false
 */
bool ModelAlgebraic::hasExtendedPattern() const
{
    return false;
}

/**
 * return an empty blockPattern if not overhead
 */
ModelAlgebraic::block_pattern_type
ModelAlgebraic::blockPattern() const
{
    return block_pattern_type( 0, 0 );
}

bool ModelAlgebraic::buildMatrixPrecond() const
{
    return !( Environment::vm()[prefixvm( this->prefix(), "preconditioner.contribution" )].as<std::string>() == "same_matrix" );
}

void ModelAlgebraic::updatePreconditioner( const vector_ptrtype& X,
                                           sparse_matrix_ptrtype& A,
                                           sparse_matrix_ptrtype& A_extended,
                                           sparse_matrix_ptrtype& Prec ) const
{
    std::string precType = option( _prefix = this->prefix(), _name = "preconditioner.contribution" ).as<std::string>();

    if ( precType == "same_matrix" )
    {
        // only copy shrared_ptr (normally already done in constructor)
        Prec = A;
    }
    else if ( precType == "standart" )
    {
        // copy standart pattern
        Prec->zero();
        Prec->addMatrix( 1., A );
    }
    else if ( precType == "extended" )
    {
        // copy standart and extended pattern
        Prec->zero();
        Prec->addMatrix( 1., A );
        if ( hasExtendedPattern() )
            Prec->addMatrix( 1., A_extended );
    }
}

ModelAlgebraic::graph_ptrtype
ModelAlgebraic::buildMatrixGraph() const
{
    return graph_ptrtype();
}

void ModelAlgebraic::updateInHousePreconditioner( sparse_matrix_ptrtype const& mat,
                                                  vector_ptrtype const& vecSol ) const
{
}

void ModelAlgebraic::updateNewtonInitialGuess( vector_ptrtype& U ) const
{
}
void ModelAlgebraic::updateJacobian( DataUpdateJacobian& data ) const
{
}
void ModelAlgebraic::updateResidual( DataUpdateResidual& data ) const
{
}
void ModelAlgebraic::updateLinearPDE( DataUpdateLinear& data ) const
{
}
void ModelAlgebraic::updatePicard( DataUpdateLinear& data ) const
{
}
double
ModelAlgebraic::updatePicardConvergence( vector_ptrtype const& Unew, vector_ptrtype const& Uold ) const
{
    return 0.;
}

} // namespace FeelModels

} // namespace Feel
