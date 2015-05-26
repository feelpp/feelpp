/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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
 \file applibase.hpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2012-01-19
 */

#ifndef APPLIBASE_HPP
#define APPLIBASE_HPP 1

#include <feel/options.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/pslogger.hpp>
#include <feel/feelcore/worldcomm.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/vectorblock.hpp>

#include <feel/feelmodels2/modelcore/feelmodelscoreconstconfig.h>
#include <feel/feelmodels2/modelcore/log.hpp>
#include <feel/feelmodels2/modelcore/timertool.hpp>


namespace Feel
{
namespace FeelModels
{


class AppliBase
{
public :
    AppliBase( std::string _theprefix,
               WorldComm const& _worldComm=Environment::worldComm(),
               std::string subPrefix="",
               std::string appliShortRepository=option(_name="exporter.directory").as<std::string>() );

    AppliBase( AppliBase const& app ) = default;
    virtual ~AppliBase();

    // worldcomm
    WorldComm const& worldComm() const;
    std::vector<WorldComm> const& worldsComm() const;
    void setWorldsComm(std::vector<WorldComm> const& _worldsComm);
    std::vector<WorldComm> const& localNonCompositeWorldsComm() const;
    void setLocalNonCompositeWorldsComm(std::vector<WorldComm> const& _worldsComm);
    virtual void createWorldsComm();
    // prefix
    std::string prefix() const;
    std::string subPrefix() const;
    // vm
    po::variables_map const& vm() const;
    // appli repository
    std::string appliRepositoryWithoutNumProc() const;
    std::string appliRepository() const;
    std::string appliShortRepository() const;
    std::string appliShortRepositoryWithNumProc() const;
    // verbose
    bool verbose() const;
    bool verboseAllProc() const;
    // info
    std::string filenameSaveInfo() const;
    void setFilenameSaveInfo(std::string s);
    virtual boost::shared_ptr<std::ostringstream> getInfo() const;
    virtual void printInfo() const;
    virtual void saveInfo() const;
    virtual void printAndSaveInfo() const;
    // timer
    TimerToolBase & timerTool(std::string s);
    void addTimerTool(std::string s,std::string fileName);
    // save assembly/solver scalability
    bool scalabilitySave() const;
    bool scalabilityReinitSaveFile() const;
    void setScalabilitySave( bool b );
    std::string scalabilityPath() const;
    void setScalabilityPath(std::string s);
    std::string scalabilityFilename() const;
    void setScalabilityFilename(std::string s);

private :
    // worldcomm
    WorldComm M_worldComm;
    std::vector<WorldComm> M_worldsComm;
    std::vector<WorldComm> M_localNonCompositeWorldsComm;
    // prefix
    std::string M_prefix;
    std::string M_subPrefix;
    // short repository name
    std::string M_appliShortRepository;
    // verbose
    bool M_verbose,M_verboseAllProc;
    // filename for save info
    std::string M_filenameSaveInfo;
    // timertool register by a name id
    std::map<std::string,std::shared_ptr<TimerToolBase> > M_mapTimerTool;
    bool M_timersActivated;
    bool M_timersSaveFileMasterRank, M_timersSaveFileMax, M_timersSaveFileMin, M_timersSaveFileMean, M_timersSaveFileAll;
    // save assembly/solver scalability
    bool M_scalabilitySave, M_scalabilityReinitSaveFile;
    std::string M_scalabilityPath;
    std::string M_scalabilityFilename;
};

//--------------------------------------------------------
// null application
struct AppliBaseNull
{
    static const bool is_class_null = true;
};

//--------------------------------------------------------
// method num application
class AppliBaseMethodsNum : public AppliBase
{
public :
    typedef AppliBase super_type;

    typedef double value_type;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_ptrtype vector_ptrtype;

    typedef backend_type::graph_ptrtype graph_ptrtype;
    typedef backend_type::indexsplit_type indexsplit_type;
    typedef backend_type::indexsplit_ptrtype indexsplit_ptrtype;

    typedef vf::BlocksBase<size_type> block_pattern_type;

    AppliBaseMethodsNum( std::string _theprefix,
                         WorldComm const& _worldComm=Environment::worldComm(),
                         std::string subPrefix="",
                         std::string appliShortRepository=option(_name="exporter.directory").as<std::string>() );

    AppliBaseMethodsNum( AppliBaseMethodsNum const& app ) = default;

    virtual ~AppliBaseMethodsNum();


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
     * return false
     */
    virtual bool hasExtendedPattern() const;

    /**
     * return an empty blockPattern if not overhead
     */
    virtual block_pattern_type blockPattern() const;

    bool buildMatrixPrecond() const;

    virtual
    void
    updatePreconditioner(const vector_ptrtype& X,
                         sparse_matrix_ptrtype& A,
                         sparse_matrix_ptrtype& A_extended,
                         sparse_matrix_ptrtype& Prec) const;

    virtual graph_ptrtype buildMatrixGraph() const;
    //----------------------------------------------------------------------------------//

    virtual void updateCLDirichlet(vector_ptrtype& U) const;// {} // const = 0;
    virtual void updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J , vector_ptrtype& R,
                                 bool BuildCstPart,
                                 sparse_matrix_ptrtype& A_extended, bool _BuildExtendedPart,
                                 bool _doClose=true, bool _doBCStrongDirichlet=true) const;// {}// = 0;
    virtual void updateResidual( const vector_ptrtype& X, vector_ptrtype& R,
                                 bool BuildCstPart, bool UseJacobianLinearTerms,
                                 bool _doClose=true, bool _doBCStrongDirichlet=true ) const;// {}// = 0;

    virtual void updateLinearPDE(const vector_ptrtype& X,sparse_matrix_ptrtype& A , vector_ptrtype& F,
                                 bool _buildCstPart,
                                 sparse_matrix_ptrtype& A_extended, bool _BuildExtendedPart,
                                 bool _doClose=true, bool _doBCStrongDirichlet=true ) const;// {}// = 0;

    //virtual void solve() {} // = 0;

    //----------------------------------------------------------------------------------//

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
};


//--------------------------------------------------------

} // namespace FeelModels
} // namespace feel


#endif //endif APPLIBASE_HPP
