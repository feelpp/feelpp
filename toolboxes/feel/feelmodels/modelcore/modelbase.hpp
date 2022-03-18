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
 \file modelbase.hpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2012-01-19
 */

#ifndef FEELPP_MODELBASE_HPP
#define FEELPP_MODELBASE_HPP 1

#include <feel/options.hpp>
#include <feel/feelcore/environment.hpp>
//#include <feel/feelcore/pslogger.hpp>
#include <feel/feelcore/worldcomm.hpp>
#include <feel/feelcore/remotedata.hpp>
#include <feel/feelmodels/modelproperties.hpp>

#include <feel/feelmodels/modelcore/feelmodelscoreconstconfig.hpp>
#include <feel/feelmodels/modelcore/log.hpp>
#include <feel/feelmodels/modelcore/timertool.hpp>

#include <feel/feelcore/tabulateinformations.hpp>

namespace Feel
{

namespace TabulateInformationTools
{
namespace FromJSON
{
tabulate_informations_ptr_t
tabulateInformationsModelFields( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );
}
}


namespace FeelModels
{

void printToolboxApplication( std::string const& toolboxName, worldcomm_t const& worldComm = Environment::worldComm() );

struct ModelBaseCommandLineOptions
{
    ModelBaseCommandLineOptions() = default;
    explicit ModelBaseCommandLineOptions( po::options_description const& _options );
    ModelBaseCommandLineOptions( ModelBaseCommandLineOptions const& ) = default;
    ModelBaseCommandLineOptions( ModelBaseCommandLineOptions && ) = default;

    po::variables_map const& vm() const
        {
            if ( M_vm.has_value() )
                return *M_vm;
            else
                return Environment::vm();
        }
private :
    std::optional<po::variables_map> M_vm;
};

struct ModelBaseRepository
{
    ModelBaseRepository( std::string const& rootDirWithoutNumProc = "", bool use_npSubDir = true, std::string const& exprRepository = "" );
    ModelBaseRepository( ModelBaseRepository const& ) = default;
    ModelBaseRepository( ModelBaseRepository && ) = default;

    std::string const& root() const { return this->rootWithNumProc(); }
    std::string const& rootWithNumProc() const { return M_rootRepositoryWithNumProc; }
    std::string const& rootWithoutNumProc() const { return M_rootRepositoryWithoutNumProc; }
    std::string const& expr() const { return M_exprRepository; }

private :
    std::string M_rootRepositoryWithNumProc, M_rootRepositoryWithoutNumProc;
    std::string M_exprRepository;
};

struct ModelBaseUpload
{
    ModelBaseUpload() = default;
    ModelBaseUpload( std::string const& desc, std::string const& basePath, worldcomm_ptr_t const& worldComm );
    ModelBaseUpload( ModelBaseUpload const& ) = default;
    ModelBaseUpload( ModelBaseUpload && ) = default;

    bool isOperational() const;

    void upload( std::string const& dataPath ) const;
    void createFolder( std::string const& folderPath, std::string const& parentId = "" ) const;

    std::string relativePath( std::string const& s ) const;

    void print() const;
private :
    void uploadPreProcess( std::string const& dataPath,
                           std::vector<std::tuple<std::string,std::time_t,std::string>> & resNewFile,
                           std::vector<std::tuple<std::string,std::time_t,std::string,std::string>> & resReplaceFile ) const;

private :
    std::shared_ptr<RemoteData> M_remoteData;
    std::string M_basePath;
    // [ folder path -> ( folder id , [ filename -> file id, last write time ] ) ]
    mutable std::map<std::string,std::pair<std::string,std::map<std::string,std::pair<std::string,std::time_t>>>> M_treeDataStructure;
};

class ModelBase : public JournalWatcher,
                  public std::enable_shared_from_this<ModelBase>
{
    using super_type = JournalWatcher;
public :
    using worldcomm_t = WorldComm;
    using worldcomm_ptr_t = std::shared_ptr<WorldComm>;
    using worldscomm_ptr_t = std::vector<std::shared_ptr<WorldComm>>;

    //!
    //! @param worldcomm communicator
    //!
    //! The worldcomm must be allocated via shared_ptr. The WorldComm can be retrieved via \c shared_from_this()
    //!
    ModelBase( std::string const& prefix, std::string const& keyword,
               worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
               std::string const& subPrefix = "",
               ModelBaseRepository const& modelRep = ModelBaseRepository(),
               ModelBaseCommandLineOptions const& modelCmdLineOpt = ModelBaseCommandLineOptions() );
    ModelBase( std::string const& prefix,
               worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
               std::string const& subPrefix = "",
               ModelBaseRepository const& modelRep = ModelBaseRepository(),
               ModelBaseCommandLineOptions const& modelCmdLineOpt = ModelBaseCommandLineOptions() )
        :
        ModelBase( prefix, prefix, worldComm, subPrefix, modelRep, modelCmdLineOpt )
        {}
    //ModelBase() : ModelBase("") {}
    ModelBase() = delete;

    ModelBase( ModelBase const& app ) = default;
    ModelBase( ModelBase && app ) = default;
    virtual ~ModelBase();

    // worldcomm
    worldcomm_ptr_t const& worldCommPtr() const;
    worldcomm_ptr_t & worldCommPtr();
    worldcomm_t & worldComm();
    worldcomm_t const& worldComm() const;
    worldscomm_ptr_t & worldsComm();
    worldscomm_ptr_t const& worldsComm() const;
    void setWorldsComm(worldscomm_ptr_t & _worldsComm);
    worldscomm_ptr_t & localNonCompositeWorldsComm() ;
    worldscomm_ptr_t const& localNonCompositeWorldsComm() const;
    void setLocalNonCompositeWorldsComm( worldscomm_ptr_t & _worldsComm);
    virtual void createWorldsComm();
    //! return variables map from command line options
    po::variables_map const& clovm() const { return M_modelCommandLineOptions.vm(); }
    // prefix
    std::string const& prefix() const;
    std::string const& subPrefix() const;
    //! keyword
    std::string const& keyword() const { return M_keyword; }
    //void setKeyword( std::string const& keyword ) { M_keyword = keyword; }
    // root repository
    ModelBaseRepository const& repository() const { return M_modelRepository; }
    std::string const& rootRepository() const;
    // verbose
    bool verbose() const;
    bool verboseAllProc() const;
    void log( std::string const& _className,std::string const& _functionName,std::string const& _msg ) const;
    // info
    void updateInformationObject( nl::json & p ) const override;

    tabulate_informations_ptr_t tabulateInformations() const;
    virtual tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const;

    virtual std::shared_ptr<std::ostringstream> getInfo() const;
    void printInfo() const { this->printInfo( this->tabulateInformations() ); }
    void saveInfo() const { this->saveInfo( this->tabulateInformations() ); }
    void printAndSaveInfo() const;
private :
    void printInfo( tabulate_informations_ptr_t const& tabInfo ) const;
    void saveInfo( tabulate_informations_ptr_t const& tabInfo ) const;
public :
    // timer
    TimerToolBase & timerTool( std::string const& s ) const;
    void addTimerTool( std::string const& s, std::string const& fileName ) const;
    // save assembly/solver scalability
    bool scalabilitySave() const;
    bool scalabilityReinitSaveFile() const;
    void setScalabilitySave( bool b );
    std::string scalabilityPath() const;
    void setScalabilityPath(std::string const& s);
    std::string scalabilityFilename() const;
    void setScalabilityFilename(std::string const& s);
    //
    bool isUpdatedForUse() const { return M_isUpdatedForUse; }
    void setIsUpdatedForUse( bool b ) { M_isUpdatedForUse = b; }
    // upload
    ModelBaseUpload const& upload() const { return M_upload; }
    void upload( std::string const& dataPath ) const;

    // model properties
    void initModelProperties();
    bool hasModelProperties() const { return (M_modelProps)? true : false; }
    std::shared_ptr<ModelProperties> modelPropertiesPtr() const { return M_modelProps; }
    ModelProperties const& modelProperties() const { return *M_modelProps; }
    ModelProperties & modelProperties() { return *M_modelProps; }
    void setModelProperties( std::shared_ptr<ModelProperties> modelProps ) { M_modelProps = modelProps; }

    /**
     * @brief Set the Model Properties object from a filename
     * 
     * @param filename file name
     */
    void setModelProperties( std::string const& filename );

    /**
     * @brief Set the Model Properties object from a json struct
     * the json may come from python 
     * @param j json data structure
     */
    void setModelProperties( nl::json const& j );

    void addParameterInModelProperties( std::string const& symbolName,double value );

    bool manageParameterValues() const { return M_manageParameterValues; }
    void setManageParameterValues( bool b ) { M_manageParameterValues = b; }
    bool manageParameterValuesOfModelProperties() const { return M_manageParameterValuesOfModelProperties; }
    void setManageParameterValuesOfModelProperties( bool b ) { M_manageParameterValuesOfModelProperties = b; }

private :
    // worldcomm
    worldcomm_ptr_t M_worldComm;
    worldscomm_ptr_t M_worldsComm;
    worldscomm_ptr_t M_localNonCompositeWorldsComm;
    // prefix
    std::string M_prefix;
    std::string M_subPrefix;
    // keyword (can be usefull in json for example)
    std::string M_keyword;
    // directory
    ModelBaseRepository M_modelRepository;
    // command line options
    ModelBaseCommandLineOptions M_modelCommandLineOptions;
    // verbose
    bool M_verbose,M_verboseAllProc;
    // timertool register by a name id
    mutable std::map<std::string,std::shared_ptr<TimerToolBase> > M_mapTimerTool;
    bool M_timersActivated;
    bool M_timersSaveFileMasterRank, M_timersSaveFileMax, M_timersSaveFileMin, M_timersSaveFileMean, M_timersSaveFileAll;
    // save assembly/solver scalability
    bool M_scalabilitySave, M_scalabilityReinitSaveFile;
    std::string M_scalabilityPath;
    std::string M_scalabilityFilename;
    //
    bool M_isUpdatedForUse;
    // upload data tools
    ModelBaseUpload M_upload;

    // model properties
    std::shared_ptr<ModelProperties> M_modelProps;
    bool M_manageParameterValues, M_manageParameterValuesOfModelProperties;

};

} // namespace FeelModels
} // namespace feel


#endif //endif FEELPP_MODELBASE_HPP
