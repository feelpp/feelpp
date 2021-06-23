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
 \file modelalgebraic.hpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2012-01-19
 */

#ifndef FEELPP_MODELALGEBRAIC_HPP
#define FEELPP_MODELALGEBRAIC_HPP 1

#include <feel/feelmodels/modelcore/modelbase.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/vectorblock.hpp>
#include <feel/feelmesh/enums.hpp>
#include <feel/feeldiscr/enums.hpp>

namespace Feel
{
namespace FeelModels
{

class ModelAlgebraic : virtual public ModelBase
{
public :
    typedef ModelBase super_type;

    typedef double value_type;
    using index_type = uint32_type;
    using size_type = index_type;
    typedef Backend<value_type,size_type> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;

    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_ptrtype vector_ptrtype;

    typedef backend_type::graph_type graph_type;
    typedef backend_type::graph_ptrtype graph_ptrtype;
    typedef backend_type::indexsplit_type indexsplit_type;
    typedef backend_type::indexsplit_ptrtype indexsplit_ptrtype;

    typedef vf::BlocksBase<size_type> block_pattern_type;

    class DataUpdateBase
    {
    public:
        DataUpdateBase() = default;
        DataUpdateBase( DataUpdateBase const& ) = default;
        DataUpdateBase( DataUpdateBase && ) = default;
        virtual ~DataUpdateBase() {}

        void addInfo( std::string const& info ) { M_infos.insert( info ); }
        void eraseInfo( std::string const& info ) { M_infos.erase( info ); }
        bool hasInfo( std::string const& info ) const { return M_infos.find( info ) != M_infos.end(); }

        void addDoubleInfo( std::string const& info, double val ) { M_doubleInfos[info] = val; }
        void eraseDoubleInfo( std::string const& info ) { M_doubleInfos.erase( info ); }
        bool hasDoubleInfo( std::string const& info ) const { return M_doubleInfos.find( info ) != M_doubleInfos.end(); }
        double doubleInfo( std::string const& info ) const
            {
                CHECK( this->hasDoubleInfo( info ) ) << "double info "<< info << "is missing";
                return M_doubleInfos.find(info)->second;
            }

        void addVectorInfo( std::string const& info, vector_ptrtype vec ) { M_vectorInfos[info] = vec; }
        void eraseVectorInfo( std::string const& info ) { M_vectorInfos.erase( info ); }
        bool hasVectorInfo( std::string const& info ) const { return M_vectorInfos.find( info ) != M_vectorInfos.end(); }
        vector_ptrtype vectorInfo( std::string const& info )
            {
                CHECK( this->hasVectorInfo( info ) ) << "vector info "<< info << "is missing";
                return M_vectorInfos.find(info)->second;
            }

        void addMatrixInfo( std::string const& info, sparse_matrix_ptrtype mat ) { M_matrixInfos[info] = mat; }
        void eraseMatrixInfo( std::string const& info ) { M_matrixInfos.erase( info ); }
        bool hasMatrixInfo( std::string const& info ) const { return M_matrixInfos.find( info ) != M_matrixInfos.end(); }
        sparse_matrix_ptrtype matrixInfo( std::string const& info )
            {
                CHECK( this->hasMatrixInfo( info ) ) << "matrix info "<< info << "is missing";
                return M_matrixInfos.find(info)->second;
            }

        void addParameterValuesInfo( std::string const& info, std::map<std::string,double> const& pv ) { M_parameterValuesInfos[info] = pv; }
        void eraseParameterValuesInfo( std::string const& info ) { M_parameterValuesInfos.erase( info ); }
        bool hasParameterValuesInfo( std::string const& info ) const { return M_parameterValuesInfos.find( info ) != M_parameterValuesInfos.end(); }
        std::map<std::string,double> & parameterValuesInfo( std::string const& info )
            {
                CHECK( this->hasParameterValuesInfo( info ) ) << "parameter values info "<< info << "is missing";
                return M_parameterValuesInfos.find(info)->second;
            }
        std::map<std::string,double> const& parameterValuesInfo( std::string const& info ) const
            {
                CHECK( this->hasParameterValuesInfo( info ) ) << "parameter values info "<< info << "is missing";
                return M_parameterValuesInfos.find(info)->second;
            }

        void copyInfos( DataUpdateBase const& dub )
            {
                for ( std::string const& info : dub.M_infos )
                    this->addInfo( info );
                for ( auto const& [info,val] : dub.M_doubleInfos )
                    this->addDoubleInfo( info,val );
                for ( auto const& [info,vec] : dub.M_vectorInfos )
                    this->addVectorInfo( info,vec );
                for ( auto const& [info,mat] : dub.M_matrixInfos )
                    this->addMatrixInfo( info,mat );
                for ( auto const& [info,pv] : dub.M_parameterValuesInfos )
                    this->addParameterValuesInfo( info,pv );
            }

        void clearInfos()
            {
                M_infos.clear();
                M_doubleInfos.clear();
                M_vectorInfos.clear();
                M_matrixInfos.clear();
                M_parameterValuesInfos.clear();
            }
    private :
        std::set<std::string> M_infos;
        std::map<std::string,double> M_doubleInfos;
        std::map<std::string,vector_ptrtype> M_vectorInfos;
        std::map<std::string,sparse_matrix_ptrtype> M_matrixInfos;
        std::map<std::string,std::map<std::string,double>> M_parameterValuesInfos;

    };

    class DataDofEliminationIds
    {
    public :
        DataDofEliminationIds() : M_hasDofEliminationIds( false ) {}
        DataDofEliminationIds( DataDofEliminationIds const& d) = default;
        DataDofEliminationIds( DataDofEliminationIds && d) = default;

        bool hasDofEliminationIds() const { return M_hasDofEliminationIds; }
        void setHasDofEliminationIds( bool b ) { M_hasDofEliminationIds = b; }
        std::set<size_type> const& dofEliminationIds() const { return M_dofEliminationIds; }
        std::set<size_type> & dofEliminationIds() { return M_dofEliminationIds; }

    private:
        bool M_hasDofEliminationIds;
        std::set<size_type> M_dofEliminationIds;
    };

    class DataUpdateLinear : public DataUpdateBase
    {
    public:
        DataUpdateLinear( const vector_ptrtype& currentSolution,
                          sparse_matrix_ptrtype matrix, vector_ptrtype rhs,
                          bool buildCstPart )
            :
            DataUpdateBase(),
            M_matrix( matrix ),
            M_rhs( rhs ),
            M_currentSolution( currentSolution ),
            M_buildCstPart( buildCstPart ),
            M_doBCStrongDirichlet( true )
            {}
        DataUpdateLinear(DataUpdateLinear const& d) = default;
        DataUpdateLinear(DataUpdateLinear && d) = default;

        sparse_matrix_ptrtype& matrix() { return M_matrix; }
        vector_ptrtype& rhs() { return M_rhs; }
        vector_ptrtype const& currentSolution() { return M_currentSolution; }
        bool buildCstPart() const { return M_buildCstPart; }
        FEELPP_DEPRECATED bool doBCStrongDirichlet() const { return M_doBCStrongDirichlet; }
        std::map<Feel::MatrixStructure,std::pair<sparse_matrix_ptrtype,double>> const& matrixToAdd() const { return M_matrixToAdd; }
        std::vector<std::pair<sparse_matrix_ptrtype,vector_ptrtype>> const& rhsToAddFromMatrixVectorProduct() const { return M_rhsToAddFromMatrixVectorProduct; }

        void setBuildCstPart( bool b ) { M_buildCstPart = b; }
        FEELPP_DEPRECATED void setDoBCStrongDirichlet( bool b ){ M_doBCStrongDirichlet = b; }
        void addMatrixToAdd( sparse_matrix_ptrtype mat, Feel::MatrixStructure matStruc, double scaling ) { M_matrixToAdd[matStruc] = std::make_pair( mat, scaling ); }
        void addRhsToAdd( sparse_matrix_ptrtype mat, vector_ptrtype vec ) { M_rhsToAddFromMatrixVectorProduct.push_back( std::make_pair(mat,vec) ); }
    private :
        sparse_matrix_ptrtype M_matrix;
        vector_ptrtype M_rhs;
        const vector_ptrtype& M_currentSolution;

        bool M_buildCstPart;
        bool M_doBCStrongDirichlet;

        std::map<Feel::MatrixStructure,std::pair<sparse_matrix_ptrtype,double>> M_matrixToAdd;
        std::vector<std::pair<sparse_matrix_ptrtype,vector_ptrtype>> M_rhsToAddFromMatrixVectorProduct;
    };

    class DataUpdateResidual : public DataUpdateBase, public DataDofEliminationIds
    {
    public:
        DataUpdateResidual( const vector_ptrtype& currentSolution, vector_ptrtype residual,
                            bool buildCstPart, bool useJacobianLinearTerms )
            :
            DataUpdateBase(),
            M_residual( residual ),
            M_currentSolution( currentSolution ),
            M_buildCstPart( buildCstPart ),
            M_useJacobianLinearTerms( useJacobianLinearTerms ),
            M_doBCStrongDirichlet( true )
            {}
        DataUpdateResidual( DataUpdateResidual const& d ) = default;
        DataUpdateResidual( DataUpdateResidual && d ) = default;

        vector_ptrtype& residual() { return M_residual; }
        vector_ptrtype const& currentSolution() { return M_currentSolution; }
        bool buildCstPart() const { return M_buildCstPart; }
        bool useJacobianLinearTerms() const { return M_useJacobianLinearTerms; }
        bool doBCStrongDirichlet() const { return M_doBCStrongDirichlet; }

        void setBuildCstPart( bool b ) { M_buildCstPart = b; }
        void setUseJacobianLinearTerms( bool b ) { M_useJacobianLinearTerms = b; }
        void setDoBCStrongDirichlet( bool b ){ M_doBCStrongDirichlet = b; }

    private :
        vector_ptrtype M_residual;
        const vector_ptrtype& M_currentSolution;
        bool M_buildCstPart;
        bool M_useJacobianLinearTerms;
        bool M_doBCStrongDirichlet;
    };

    class DataUpdateJacobian : public DataUpdateBase, public DataDofEliminationIds
    {
    public:
        DataUpdateJacobian( const vector_ptrtype& currentSolution, sparse_matrix_ptrtype jacobian,
                            vector_ptrtype vectorUsedInStrongDirichlet, bool buildCstPart )
            :
            DataUpdateBase(),
            M_jacobian( jacobian ),
            M_vectorUsedInStrongDirichlet( vectorUsedInStrongDirichlet ),
            M_currentSolution( currentSolution ),
            M_buildCstPart( buildCstPart ),
            M_doBCStrongDirichlet( true )
            {}

        DataUpdateJacobian( DataUpdateJacobian const& d) = default;
        DataUpdateJacobian( DataUpdateJacobian && d) = default;

        sparse_matrix_ptrtype& jacobian() { return M_jacobian; }
        vector_ptrtype& vectorUsedInStrongDirichlet() { return M_vectorUsedInStrongDirichlet; }
        vector_ptrtype const& currentSolution() { return M_currentSolution; }

        bool buildCstPart() const { return M_buildCstPart; }
        FEELPP_DEPRECATED bool doBCStrongDirichlet() const { return M_doBCStrongDirichlet; }

        void setBuildCstPart( bool b ) { M_buildCstPart = b; }
        void setDoBCStrongDirichlet( bool b ){ M_doBCStrongDirichlet = b; }

    private :
        sparse_matrix_ptrtype M_jacobian;
        vector_ptrtype M_vectorUsedInStrongDirichlet;
        const vector_ptrtype& M_currentSolution;
        bool M_buildCstPart;
        bool M_doBCStrongDirichlet;
    };

    class DataDofEliminationIdsByEntity
    {
    public :
        DataDofEliminationIdsByEntity() = default;
        DataDofEliminationIdsByEntity( DataDofEliminationIdsByEntity const& d) = default;
        DataDofEliminationIdsByEntity( DataDofEliminationIdsByEntity && d) = default;

        bool hasDofEliminationIds( ElementsType e ) const { return M_dofEliminationIds.find( e ) != M_dofEliminationIds.end(); }
        bool hasDofEliminationIds() const
            {
                return this->hasDofEliminationIds( MESH_ELEMENTS ) ||
                    this->hasDofEliminationIds( MESH_FACES ) ||
                    this->hasDofEliminationIds( MESH_EDGES ) ||
                    this->hasDofEliminationIds( MESH_POINTS );
            }
        std::set<size_type> const& dofEliminationIds( ElementsType e ) const
            {
                CHECK( this->hasDofEliminationIds(e) ) << "entity not registered";
                return M_dofEliminationIds.find( e )->second;
            }
        std::set<size_type> & dofEliminationIds( ElementsType e ) { return M_dofEliminationIds[e]; }
        void initDofEliminationIds( ElementsType e ) { M_dofEliminationIds[e]; }
        void addDofEliminationIds( ElementsType e, size_type id ) { M_dofEliminationIds[e].insert( id ); }
        void addDofEliminationIds( ElementsType e, std::set<size_type> const& ids ) { M_dofEliminationIds[e].insert( ids.begin(), ids.end() ); }
    private:
        std::map<ElementsType,std::set<size_type>> M_dofEliminationIds;
    };

    class DataNewtonInitialGuess : public DataUpdateBase, public DataDofEliminationIdsByEntity
    {
    public:
        DataNewtonInitialGuess( vector_ptrtype& initialGuess )
            :
            DataUpdateBase(),
            DataDofEliminationIdsByEntity(),
            M_initialGuess( initialGuess )
            {}
        DataNewtonInitialGuess( DataNewtonInitialGuess const& d) = default;
        DataNewtonInitialGuess( DataNewtonInitialGuess && d) = default;
        vector_ptrtype & initialGuess() { return M_initialGuess; }


    private:
        vector_ptrtype& M_initialGuess;
    };

    ModelAlgebraic( std::string _theprefix, std::string const& keyword,
                    worldcomm_ptr_t const& _worldComm=Environment::worldCommPtr(),
                    std::string const& subPrefix="",
                    ModelBaseRepository const& modelRep = ModelBaseRepository(),
                    ModelBaseCommandLineOptions const& modelCmdLineOpt = ModelBaseCommandLineOptions() );
    ModelAlgebraic( std::string _theprefix,
                    worldcomm_ptr_t const& _worldComm=Environment::worldCommPtr(),
                    std::string const& subPrefix="",
                    ModelBaseRepository const& modelRep = ModelBaseRepository(),
                    ModelBaseCommandLineOptions const& modelCmdLineOpt = ModelBaseCommandLineOptions() )
        :
        ModelAlgebraic( _theprefix, _theprefix, _worldComm, subPrefix, modelRep, modelCmdLineOpt )
        {}

    ModelAlgebraic( ModelAlgebraic const& app ) = default;

    virtual ~ModelAlgebraic();


    // verbose
    bool verboseSolverTimer() const;
    bool verboseSolverTimerAllProc() const;
    // do rebuild cst part in linear/jacobian or use jac for residual
    bool rebuildCstPartInLinearSystem() const;
    void setRebuildCstPartInLinearSystem(bool b);
    bool useLinearJacobianInResidual() const;
    void setUseLinearJacobianInResidual(bool b);
    bool rebuildLinearPartInJacobian() const;
    void setRebuildLinearPartInJacobian(bool b);
    // a utiliser avec precaution!!!
    bool rebuildCstPartInResidual() const;
    void setRebuildCstPartInResidual(bool b);
    // define an other matrix/vector to store the cst part
    bool useCstMatrix() const;
    void setUseCstMatrix(bool b);
    bool useCstVector() const;
    void setUseCstVector(bool b);
    // allow to rebuild cst part (once at next solve) if some parameters (model,time mode,..) change
    bool needToRebuildCstPart() const;
    void setNeedToRebuildCstPart(bool b);
    // an option
    bool errorIfSolverNotConverged() const;
    void setErrorIfSolverNotConverged( bool b );
    // save a python script to view graph
    bool printGraph() const;
    void setPrintGraph(bool b);
    std::string printGraphFileName() const;
    void setPrintGraphFileName(std::string s);

    //----------------------------------------------------------------------------------//

    /**
     * return an empty blockPattern if not overhead
     */
    virtual block_pattern_type blockPattern() const;

    virtual void updateInHousePreconditioner( DataUpdateLinear & data ) const;
    virtual void updateInHousePreconditioner( DataUpdateJacobian & data ) const;

    virtual BlocksBaseGraphCSR buildBlockMatrixGraph() const;
    virtual graph_ptrtype buildMatrixGraph() const;

    //----------------------------------------------------------------------------------//

    virtual void updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const;
    virtual void updateJacobian( DataUpdateJacobian & data ) const;
    virtual void updateJacobianDofElimination( DataUpdateJacobian & data ) const;
    virtual void updateResidual( DataUpdateResidual & data ) const;
    virtual void updateResidualDofElimination( DataUpdateResidual & data ) const;
    virtual void updateLinearPDE( DataUpdateLinear & data ) const;
    virtual void updateLinearPDEDofElimination( DataUpdateLinear & data ) const;

    //----------------------------------------------------------------------------------//
    virtual void preSolveNewton( vector_ptrtype rhs, vector_ptrtype sol ) const {}
    virtual void postSolveNewton( vector_ptrtype rhs, vector_ptrtype sol ) const {}
    virtual void preSolvePicard( vector_ptrtype rhs, vector_ptrtype sol ) const {}
    virtual void postSolvePicard( vector_ptrtype rhs, vector_ptrtype sol ) const {}
    virtual void preSolveLinear( vector_ptrtype rhs, vector_ptrtype sol ) const {}
    virtual void postSolveLinear( vector_ptrtype rhs, vector_ptrtype sol ) const {}
    //----------------------------------------------------------------------------------//
    virtual void updateNewtonIteration( int step, vector_ptrtype residual, vector_ptrtype sol, typename backend_type::solvernonlinear_type::UpdateIterationData const& data ) const {}
    virtual void updatePicardIteration( int step, vector_ptrtype sol ) const {}

    //! index start of (sub-)block
    size_type rowStartInMatrix() const { return this->startBlockSpaceIndexMatrixRow(); }
    size_type colStartInMatrix() const { return this->startBlockSpaceIndexMatrixCol(); }
    size_type rowStartInVector() const { return this->startBlockSpaceIndexVector(); }
    size_type startBlockSpaceIndexMatrixRow() const { return M_startBlockSpaceIndexMatrixRow; }
    size_type startBlockSpaceIndexMatrixCol() const { return M_startBlockSpaceIndexMatrixCol; }
    size_type startBlockSpaceIndexVector() const { return M_startBlockSpaceIndexVector; }
    void setStartBlockSpaceIndexMatrixRow( size_type s ) { M_startBlockSpaceIndexMatrixRow = s; }
    void setStartBlockSpaceIndexMatrixCol( size_type s ) { M_startBlockSpaceIndexMatrixCol = s; }
    void setStartBlockSpaceIndexVector( size_type s ) { M_startBlockSpaceIndexVector = s; }
    void setStartBlockSpaceIndex( size_type s ) { this->setStartBlockSpaceIndexMatrixRow( s ); this->setStartBlockSpaceIndexMatrixCol( s ); this->setStartBlockSpaceIndexVector( s ); }

    size_type startSubBlockSpaceIndex( std::string const& name ) const
        {
            auto itFind = M_startSubBlockSpaceIndex.find( name );
            if ( itFind != M_startSubBlockSpaceIndex.end() )
                return itFind->second;
            return invalid_v<size_type>;
        }
    bool hasStartSubBlockSpaceIndex( std::string const& name ) const { return (this->startSubBlockSpaceIndex( name ) != invalid_v<size_type>); }
    void setStartSubBlockSpaceIndex( std::string const& name, size_type s ) { M_startSubBlockSpaceIndex[name] = s; }

    //! update data usefull for mpi synchronization of NewtonInitialGuess, impose value in residual or jacobian
    template <typename DataType>
    void updateDofEliminationIds( std::string const& spaceName, DataType & data ) const
        {
            auto itFindDofEliminationIds = M_dofEliminationIds.find( spaceName );
            if ( itFindDofEliminationIds != M_dofEliminationIds.end() )
            {
                auto const& dofIds = itFindDofEliminationIds->second;
                this->updateDofEliminationIds( spaceName, dofIds, data );
            }
        }
    void updateDofEliminationIds( std::string const& spaceName, std::map<ElementsType, std::tuple<std::set<size_type>,std::set<size_type>>> const& dofIds, DataNewtonInitialGuess & data ) const;
    void updateDofEliminationIds( std::string const& spaceName, std::map<ElementsType, std::tuple<std::set<size_type>,std::set<size_type>>> const& dofIds, DataUpdateResidual & data ) const;
    void updateDofEliminationIds( std::string const& spaceName, std::map<ElementsType, std::tuple<std::set<size_type>,std::set<size_type>>> const& dofIds, DataUpdateJacobian & data ) const;
    std::set<size_type> & dofEliminationIdsAll( std::string const& spaceName, ElementsType e ) { return std::get<0>( M_dofEliminationIds[spaceName][e] ); }
    std::set<size_type> & dofEliminationIdsMultiProcess( std::string const& spaceName, ElementsType e ) { return std::get<1>( M_dofEliminationIds[spaceName][e] ); }

    std::map<std::string,std::map<ElementsType, std::tuple<std::set<size_type>,std::set<size_type>>>> const& dofEliminationIds() const { return M_dofEliminationIds; }
    std::map<ElementsType, std::tuple<std::set<size_type>,std::set<size_type>>> const& dofEliminationIds( std::string const& spaceName ) const
        {
            CHECK( this->hasDofEliminationIds( spaceName ) ) << "no space name registered : " << spaceName;
            return M_dofEliminationIds.find( spaceName )->second;
        }
    bool hasDofEliminationIds( std::string const& spaceName ) const { return M_dofEliminationIds.find( spaceName ) != M_dofEliminationIds.end(); }

    template <typename SpaceType, typename RangeType>
    void updateDofEliminationIds( std::string const& spaceName, std::shared_ptr<SpaceType> thespace, RangeType const& therange, ComponentType c1 = ComponentType::NO_COMPONENT )
        {
            ElementsType et = (ElementsType)boost::get<0>( therange ).value;
            auto dofsToAdd = thespace->dofs( therange, c1 );
            thespace->dof()->updateIndexSetWithParallelMissingDof( dofsToAdd );
            this->dofEliminationIdsAll(spaceName,et).insert( dofsToAdd.begin(), dofsToAdd.end() );
            auto dofsMultiProcessToAdd = thespace->dofs( therange, c1, true );
            this->dofEliminationIdsMultiProcess(spaceName,et).insert( dofsMultiProcessToAdd.begin(), dofsMultiProcessToAdd.end() );
        }
private :
    // verbose
    bool M_verboseSolverTimer,M_verboseSolverTimerAllProc;

    bool M_rebuildCstPartInLinearSystem;
    bool M_useLinearJacobianInResidual;
    bool M_rebuildLinearPartInJacobian;
    bool M_rebuildCstPartInResidual;
    bool M_useCstMatrix,M_useCstVector;
    bool M_needToRebuildCstPart;

    bool M_errorIfSolverNotConverged;
    // save a python script to view graph
    bool M_printGraph;
    std::string M_printGraphFileName;

    //! index start of (sub-)block
    size_type M_startBlockSpaceIndexMatrixRow, M_startBlockSpaceIndexMatrixCol, M_startBlockSpaceIndexVector;
    std::map<std::string,size_type> M_startSubBlockSpaceIndex;
    //! dofs eliminiation ( spaceName -> ( ElementsType -> ( all dofs, only dofs at interprocess that the value can be used) ) )
    std::map<std::string,std::map<ElementsType, std::tuple<std::set<size_type>,std::set<size_type> > > > M_dofEliminationIds;


};


} // namespace FeelModels
} // namespace feel


#endif // FEELPP_MODELALGEBRAIC_HPP
