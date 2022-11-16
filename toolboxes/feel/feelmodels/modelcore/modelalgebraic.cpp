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

namespace Feel {

namespace FeelModels {


ModelAlgebraic::ModelAlgebraic( std::string _theprefix, std::string const& keyword,
                                worldcomm_ptr_t const& _worldComm,
                                std::string const& subPrefix,
                                ModelBaseRepository const& modelRep,
                                ModelBaseCommandLineOptions const& modelCmdLineOpt )
    :
    super_type( _theprefix,keyword,_worldComm,subPrefix,modelRep, modelCmdLineOpt ),
    M_verboseSolverTimer( boption(_name="verbose_solvertimer",_prefix=this->prefix(),_vm=this->clovm()) ),
    M_verboseSolverTimerAllProc( boption(_name="verbose_solvertimer_allproc",_prefix=this->prefix(),_vm=this->clovm()) ),
    M_rebuildCstPartInLinearSystem( boption(_name="linearsystem-cst-update",_prefix=this->prefix(),_vm=this->clovm()) ),
    M_useLinearJacobianInResidual( boption(_name="residual-uselinearjac",_prefix=this->prefix(),_vm=this->clovm()) ),
    M_rebuildLinearPartInJacobian( boption(_name="jacobian-linear-update",_prefix=this->prefix(),_vm=this->clovm()) ),
    M_rebuildCstPartInResidual(true), // not an option (just opitmisation with semi-implicit)
    M_useCstMatrix( boption(_name="use-cst-matrix",_prefix=this->prefix(),_vm=this->clovm()) ),
    M_useCstVector( boption(_name="use-cst-vector",_prefix=this->prefix(),_vm=this->clovm()) ),
    M_needToRebuildCstPart( false ),
    M_errorIfSolverNotConverged( boption(_name="error-if-solver-not-converged",_prefix=this->prefix(),_vm=this->clovm()) ),
    M_printGraph( boption(_name="graph-print-python",_prefix=this->prefix(),_vm=this->clovm()) ),
    M_startBlockSpaceIndexMatrixRow(0),
    M_startBlockSpaceIndexMatrixCol(0),
    M_startBlockSpaceIndexVector(0)
{

    //-----------------------------------------------------------------------//
    //-----------------------------------------------------------------------//
    if (this->clovm().count(prefixvm(this->prefix(),"graph-print-python-filename")))
        M_printGraphFileName = soption(_name="graph-print-python-filename",_prefix=this->prefix(),_vm=this->clovm());
    else
        M_printGraphFileName = this->prefix()+".graphPython.py";

}
ModelAlgebraic::~ModelAlgebraic()
{}



// verbose
bool
ModelAlgebraic::verboseSolverTimer() const { return M_verboseSolverTimer; }
bool
ModelAlgebraic::verboseSolverTimerAllProc() const { return M_verboseSolverTimerAllProc; }
// do rebuild cst part in linear/jacobian or use jac for residual
bool
ModelAlgebraic::rebuildCstPartInLinearSystem() const { return M_rebuildCstPartInLinearSystem; }
void
ModelAlgebraic::setRebuildCstPartInLinearSystem(bool b) { M_rebuildCstPartInLinearSystem=b; }
bool
ModelAlgebraic::useLinearJacobianInResidual() const { return M_useLinearJacobianInResidual; }
void
ModelAlgebraic::setUseLinearJacobianInResidual(bool b) { M_useLinearJacobianInResidual=b; }
bool
ModelAlgebraic::rebuildLinearPartInJacobian() const { return M_rebuildLinearPartInJacobian; }
void
ModelAlgebraic::setRebuildLinearPartInJacobian(bool b) { M_rebuildLinearPartInJacobian=b; }
// a utiliser avec precaution!!!
bool
ModelAlgebraic::rebuildCstPartInResidual() const { return M_rebuildCstPartInResidual; }
void
ModelAlgebraic::setRebuildCstPartInResidual(bool b) { M_rebuildCstPartInResidual=b; }
// define an other matrix/vector to store the cst part
bool
ModelAlgebraic::useCstMatrix() const { return M_useCstMatrix; }
void
ModelAlgebraic::setUseCstMatrix(bool b) { M_useCstMatrix = b; }
bool
ModelAlgebraic::useCstVector() const { return M_useCstVector; }
void
ModelAlgebraic::setUseCstVector(bool b) { M_useCstVector = b; }
// allow to rebuild cst part (once at next solve) if some parameters (model,time mode,..) change
bool
ModelAlgebraic::needToRebuildCstPart() const { return M_needToRebuildCstPart; }
void
ModelAlgebraic::setNeedToRebuildCstPart(bool b) { M_needToRebuildCstPart = b; }
// an option
bool
ModelAlgebraic::errorIfSolverNotConverged() const { return M_errorIfSolverNotConverged; }
void
ModelAlgebraic::setErrorIfSolverNotConverged( bool b ) { M_errorIfSolverNotConverged= b; }
// save a python script to view graph
bool
ModelAlgebraic::printGraph() const { return M_printGraph; }
void
ModelAlgebraic::setPrintGraph(bool b) { M_printGraph=b; }
std::string
ModelAlgebraic::printGraphFileName() const { return M_printGraphFileName; }
void
ModelAlgebraic::setPrintGraphFileName(std::string s) { M_printGraphFileName=s; }

/**
 * return an empty blockPattern if not overhead
 */
ModelAlgebraic::block_pattern_type
ModelAlgebraic::blockPattern() const
{
    return block_pattern_type(0,0);
}

BlocksBaseGraphCSR
ModelAlgebraic::buildBlockMatrixGraph() const
{
    BlocksBaseGraphCSR myblockGraph(0,0);
    return myblockGraph;

}

ModelAlgebraic::graph_ptrtype
ModelAlgebraic::buildMatrixGraph() const
{
    auto blockGraph = this->buildBlockMatrixGraph();
    if ( blockGraph.nRow() == 0 || blockGraph.nCol() == 0 )
        return graph_ptrtype();

    if ( blockGraph.nRow() == 1 && blockGraph.nCol() == 1 )
        return blockGraph(0,0);
    else
        return graph_ptrtype( new graph_type( blockGraph ) );
}

void
ModelAlgebraic::updateInHousePreconditioner( DataUpdateLinear & data ) const
{}
void
ModelAlgebraic::updateInHousePreconditioner( DataUpdateJacobian & data ) const
{}

void
ModelAlgebraic::updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const
{}
void
ModelAlgebraic::updateJacobian( DataUpdateJacobian & data ) const
{}
void
ModelAlgebraic::updateJacobianDofElimination( DataUpdateJacobian & data ) const
{}
void
ModelAlgebraic::updateResidual( DataUpdateResidual & data ) const
{}
void
ModelAlgebraic::updateResidualDofElimination( DataUpdateResidual & data ) const
{}
void
ModelAlgebraic::updateLinearPDE( DataUpdateLinear & data ) const
{}
void
ModelAlgebraic::updateLinearPDEDofElimination( DataUpdateLinear & data ) const
{}


void
ModelAlgebraic::updateDofEliminationIds( std::string const& spaceName,
                                         std::map<ElementsType, std::tuple<std::set<size_type>,std::set<size_type>>> const& dofIds,
                                         DataNewtonInitialGuess & data ) const
{
    CHECK( this->hasStartSubBlockSpaceIndex( spaceName ) ) << "no space name registered : " << spaceName;
    int spaceIndexVector = this->startBlockSpaceIndexVector() + this->startSubBlockSpaceIndex( spaceName );
    std::vector<ElementsType> fromEntities = { MESH_ELEMENTS, MESH_FACES, MESH_EDGES, MESH_POINTS };
    auto dm = data.initialGuess()->mapPtr();
    for ( ElementsType entity : fromEntities )
    {
        auto itFindDofIdsByEntity = dofIds.find( entity );
        if ( itFindDofIdsByEntity != dofIds.end() )
            dm->dofIdToContainerId( spaceIndexVector, std::get<0>( itFindDofIdsByEntity->second ), //std::get<1>( itFindDofIdsByEntity->second ),
                                    data.dofEliminationIds( entity ) );
    }
}
void
ModelAlgebraic::updateDofEliminationIds( std::string const& spaceName,
                                         std::map<ElementsType, std::tuple<std::set<size_type>,std::set<size_type>>> const& dofIds,
                                         DataUpdateResidual & data ) const
{
    CHECK( this->hasStartSubBlockSpaceIndex( spaceName ) ) << "no space name registered : " << spaceName;
    int spaceIndexVector = this->startBlockSpaceIndexVector() + this->startSubBlockSpaceIndex( spaceName );
    std::vector<ElementsType> fromEntities = { MESH_ELEMENTS, MESH_FACES, MESH_EDGES, MESH_POINTS };
    auto dm = data.residual()->mapPtr();
    for ( ElementsType entity : fromEntities )
    {
        auto itFindDofIdsByEntity = dofIds.find( entity );
        if ( itFindDofIdsByEntity != dofIds.end() )
        {
            data.setHasDofEliminationIds( true );
            dm->dofIdToContainerId( spaceIndexVector,std::get<0>( itFindDofIdsByEntity->second ),
                                    data.dofEliminationIds() );
        }
    }
}
void
ModelAlgebraic::updateDofEliminationIds( std::string const& spaceName,
                                         std::map<ElementsType, std::tuple<std::set<size_type>,std::set<size_type>>> const& dofIds,
                                         DataUpdateJacobian & data ) const
{
    CHECK( this->hasStartSubBlockSpaceIndex( spaceName ) ) << "no space name registered : " << spaceName;
    int spaceIndexVector = this->startBlockSpaceIndexVector() + this->startSubBlockSpaceIndex( spaceName );
    std::vector<ElementsType> fromEntities = { MESH_ELEMENTS, MESH_FACES, MESH_EDGES, MESH_POINTS };
    auto dm = data.jacobian()->mapRowPtr();
    for ( ElementsType entity : fromEntities )
    {
        auto itFindDofIdsByEntity = dofIds.find( entity );
        if ( itFindDofIdsByEntity != dofIds.end() )
        {
            data.setHasDofEliminationIds( true );
            dm->dofIdToContainerId( spaceIndexVector,std::get<0>( itFindDofIdsByEntity->second ),
                                    data.dofEliminationIds() );
        }
    }
}



} // namespace FeelModels

} // namespace Feel
