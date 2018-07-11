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
 \file modelbase.cpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2012-01-19
 */

#include <feel/feelmodels/modelcore/modelbase.hpp>

namespace Feel {

namespace FeelModels {

namespace ToolboxesDetail
{
void removeTrailingSlash( std::string & s )
{
    if ( Feel::filename_is_dot( fs::path( s ).filename() ) )
        s = fs::path( s ).parent_path().string();
}
}

ModelBaseRepository::ModelBaseRepository( std::string const& rootDirWithoutNumProc )
{
    if ( rootDirWithoutNumProc.empty() )
    {
        M_rootRepositoryWithoutNumProc = Environment::appRepositoryWithoutNumProc();
        M_rootRepositoryWithNumProc = Environment::appRepository();
        M_exprRepository = Environment::exprRepository();
    }
    else
    {
        M_rootRepositoryWithoutNumProc = Environment::expand( rootDirWithoutNumProc );
        if ( fs::path( M_rootRepositoryWithoutNumProc ).is_relative() )
            M_rootRepositoryWithoutNumProc = (fs::path(Environment::rootRepository())/fs::path(M_rootRepositoryWithoutNumProc)).string();
        std::string npSubDir = (boost::format( "np_%1%" ) % Environment::worldComm().localSize() ).str();
        M_rootRepositoryWithNumProc = ( fs::path(M_rootRepositoryWithoutNumProc) / fs::path( npSubDir ) ).string();
        M_exprRepository = ( fs::path(M_rootRepositoryWithoutNumProc) / fs::path( "exprs" ) ).string();
    }

    ToolboxesDetail::removeTrailingSlash( M_rootRepositoryWithoutNumProc );
    ToolboxesDetail::removeTrailingSlash( M_rootRepositoryWithNumProc );
    ToolboxesDetail::removeTrailingSlash( M_exprRepository );
}

ModelBaseUpload::ModelBaseUpload( std::string const& desc, std::string const& basePath, WorldComm const& worldComm )
    :
    M_basePath( basePath )
{
    if ( desc.empty() )
        return;
    M_remoteData.reset( new RemoteData( desc,worldComm ) );
}

bool
ModelBaseUpload::isOperational() const
{
    if ( !M_remoteData )
        return false;
    return M_remoteData->canUpload();
}



void
ModelBaseUpload::uploadPreProcess( std::string const& dataPath,
                                   std::vector<std::tuple<std::string,std::time_t,std::string>> & resNewFile,
                                   std::vector<std::tuple<std::string,std::time_t,std::string,std::string>> & resReplaceFile ) const
{
    fs::path dataFsPath( dataPath );
    // if dir, loop on all files/dir inside
    if ( fs::is_directory( dataFsPath ) )
    {
        fs::recursive_directory_iterator end_itr;
        for (fs::recursive_directory_iterator itr( dataFsPath ); itr != end_itr; ++itr)
        {
            fs::path fileFsPath = itr->path();
            if ( !fs::is_regular_file( fileFsPath ) )
                continue;
            if ( fs::is_symlink( fileFsPath ) ) // ignore symlink
                continue;
            uploadPreProcess( fileFsPath.string(), resNewFile, resReplaceFile );
        }
        return;
    }

    std::string relDataPath = this->relativePath( dataPath );
    CHECK( !relDataPath.empty() ) << "not relative path";

    fs::path folder = fs::path( relDataPath ).parent_path();

    std::time_t dataLWT = fs::last_write_time( dataFsPath );

    std::string parentId;
    auto itFindFolder = M_treeDataStructure.find( folder.string() );
    if ( itFindFolder != M_treeDataStructure.end() )
    {
        parentId = itFindFolder->second.first;
        auto const& filesInFolder = itFindFolder->second.second;
        auto itFindFile = filesInFolder.find( dataFsPath.filename().string() );
        if ( itFindFile != filesInFolder.end() )
        {
            auto const& fileInfo = itFindFile->second;
            if ( fileInfo.second == fs::last_write_time( dataFsPath ) )
                return;
            else
                resReplaceFile.push_back( std::make_tuple(dataPath, dataLWT, folder.string(), fileInfo.first) );
        }
        else
        {
            resNewFile.push_back( std::make_tuple(dataPath, dataLWT, folder.string()) );
        }
    }
    else
    {
        auto newFolders = M_remoteData->createFolder( folder.string(), "", false );
        std::string curFolderInTree;
        for ( int k=0;k<newFolders.size();++k )
        {
            std::string curFolderName = newFolders[k].first;
            std::string curFolderId = newFolders[k].second;
            parentId = curFolderId;
            curFolderInTree = (k==0)? curFolderName : (fs::path(curFolderInTree)/curFolderName).string();
            auto itFindFolder2 = M_treeDataStructure.find( curFolderInTree );
            if ( itFindFolder2 == M_treeDataStructure.end() )
                M_treeDataStructure[ curFolderInTree ] = std::make_pair( curFolderId,std::map<std::string,std::pair<std::string,std::time_t>>() );
        }
        resNewFile.push_back( std::make_tuple(dataPath, dataLWT, folder.string()/*parentId*/) );
    }


}

void
ModelBaseUpload::upload( std::string const& dataPath ) const
{
    if ( !this->isOperational() )
        return;

    if ( M_remoteData->worldComm().isMasterRank() )
    {
        // prepare datas to upload and create folders on remote server
        std::vector<std::tuple<std::string,std::time_t,std::string>> dataPreProcessNewFile;
        std::vector<std::tuple<std::string,std::time_t,std::string,std::string>> dataPreProcessReplaceFile;
        this->uploadPreProcess( dataPath, dataPreProcessNewFile, dataPreProcessReplaceFile );

        if ( dataPreProcessNewFile.empty() && dataPreProcessReplaceFile.empty() )
            return;
        // generate data to upload container
        std::vector<std::pair<std::string,std::string> > dataToUpload( dataPreProcessNewFile.size() );
        for ( int k=0;k<dataPreProcessNewFile.size();++k )
        {
            std::string const& folder = std::get<2>( dataPreProcessNewFile[k] );
            auto itFindFolder = M_treeDataStructure.find( folder );
            CHECK( itFindFolder != M_treeDataStructure.end() ) << "folder not found";
            std::string const& parentId = itFindFolder->second.first;
            dataToUpload[k] = std::make_pair( std::get<0>( dataPreProcessNewFile[k] ), parentId );
        }
        // generate files to replace container
        std::vector<std::pair<std::string,std::string> > filesToReplace( dataPreProcessReplaceFile.size() );
        for ( int k=0;k<dataPreProcessReplaceFile.size();++k )
        {
            std::string const& filePath = std::get<0>( dataPreProcessReplaceFile[k] );
            std::string const& fileId = std::get<3>( dataPreProcessReplaceFile[k] );
            filesToReplace[k] = std::make_pair( filePath,fileId );
        }
        // apply data upload on remote server
        auto resUpload = M_remoteData->upload( dataToUpload, false );
        CHECK( resUpload.size() == dataToUpload.size() ) << "failure in upload";
        // update tree data structure
        for ( int k=0;k<dataPreProcessNewFile.size();++k )
        {
            CHECK( resUpload[k].size() == 1 ) << "must be 1";
            std::string const& fileIdUploaded = resUpload[k][0];
            std::string const& dataPathUploaded = std::get<0>( dataPreProcessNewFile[k] );
            std::time_t lwt = std::get<1>( dataPreProcessNewFile[k] );
            std::string const& folder = std::get<2>( dataPreProcessNewFile[k] );

            std::string relDataPathUploaded = this->relativePath( dataPathUploaded );
            CHECK( !relDataPathUploaded.empty() ) << "not relative path";
            M_treeDataStructure[folder].second[fs::path( relDataPathUploaded ).filename().string()] = std::make_pair(fileIdUploaded,lwt);
        }

        // apply replace files remote server
        M_remoteData->replaceFile( filesToReplace );
        // update tree data structure
        for ( int k=0;k<dataPreProcessReplaceFile.size();++k )
        {
            std::string const& dataPathUploaded = std::get<0>( dataPreProcessReplaceFile[k] );
            std::time_t lwt = std::get<1>( dataPreProcessReplaceFile[k] );
            std::string const& folder = std::get<2>( dataPreProcessReplaceFile[k] );

            std::string relDataPathUploaded = this->relativePath( dataPathUploaded );
            CHECK( !relDataPathUploaded.empty() ) << "not relative path";
            M_treeDataStructure[folder].second[fs::path( relDataPathUploaded ).filename().string()].second = lwt;
        }

        //this->print();
    }
}
void
ModelBaseUpload::createFolder( std::string const& folderPath, std::string const& parentId ) const
{
    if ( this->isOperational() )
        M_remoteData->createFolder( folderPath, parentId, false );
}
std::string
ModelBaseUpload::relativePath( std::string const& s ) const
{
    fs::path relPath = fs::relative( s,M_basePath );
    return relPath.string();
}

void
ModelBaseUpload::print() const
{
    std::cout << "ModelBaseUpload::print\n";
    for ( auto const& folder : M_treeDataStructure )
    {
        auto const& folderData = folder.second;
        std::cout << folder.first << " [id=" << folderData.first << "]\n";
        auto const& files = folderData.second;
        for ( auto const& fileData : files )
            std::cout << "  -- " << fileData.first << " ["<< fileData.second << "]\n";
    }
}


ModelBase::ModelBase( std::string const& prefix,
                      WorldComm const& worldComm,
                      std::string const& subPrefix,
                      ModelBaseRepository const& modelRep )
    :
    M_worldComm(worldComm),
    M_worldsComm( std::vector<WorldComm>(1,worldComm) ),
    M_localNonCompositeWorldsComm( std::vector<WorldComm>(1,worldComm) ),
    M_prefix( prefix ),
    M_subPrefix( subPrefix ),
    M_modelRepository( modelRep ),
    M_verbose( boption(_name="verbose",_prefix=this->prefix()) ),
    M_verboseAllProc( boption(_name="verbose_allproc",_prefix=this->prefix()) ),
    M_filenameSaveInfo( prefixvm(this->prefix(),prefixvm(this->subPrefix(),"toolbox-info.txt")) ),
    M_timersActivated( boption(_name="timers.activated",_prefix=this->prefix()) ),
    M_timersSaveFileMasterRank( boption(_name="timers.save-master-rank",_prefix=this->prefix()) ),
    M_timersSaveFileMax( boption(_name="timers.save-max",_prefix=this->prefix()) ),
    M_timersSaveFileMin( boption(_name="timers.save-min",_prefix=this->prefix()) ),
    M_timersSaveFileMean( boption(_name="timers.save-mean",_prefix=this->prefix()) ),
    M_timersSaveFileAll( boption(_name="timers.save-all",_prefix=this->prefix()) ),
    M_scalabilitySave( boption(_name="scalability-save",_prefix=this->prefix()) ),
    M_scalabilityReinitSaveFile( boption(_name="scalability-reinit-savefile",_prefix=this->prefix()) ),
    M_upload( soption(_name="upload",_prefix=this->prefix()), this->repository().rootWithoutNumProc(), M_worldComm )
{
    if (Environment::vm().count(prefixvm(this->prefix(),"scalability-path")))
        M_scalabilityPath = Environment::vm()[prefixvm(this->prefix(),"scalability-path")].as< std::string >();
    else
        M_scalabilityPath = this->repository().rootWithoutNumProc();

    if (Environment::vm().count(prefixvm(this->prefix(),"scalability-filename")))
        M_scalabilityFilename = Environment::vm()[prefixvm(this->prefix(),"scalability-filename")].as< std::string >();
    else
        M_scalabilityFilename = this->prefix()+".scalibility";

    if ( M_upload.isOperational() )
    {
        fs::path dirRelative = fs::relative( this->repository().rootWithNumProc(),
                                             this->repository().rootWithoutNumProc() );
        M_upload.createFolder( dirRelative.string() );
    }

}
ModelBase::~ModelBase()
{}

WorldComm const&
ModelBase::worldComm() const { return M_worldComm; }
std::vector<WorldComm> const&
ModelBase::worldsComm() const { return M_worldsComm; }
void
ModelBase::setWorldsComm(std::vector<WorldComm> const& _worldsComm) { M_worldsComm=_worldsComm; }
std::vector<WorldComm> const&
ModelBase::localNonCompositeWorldsComm() const { return M_localNonCompositeWorldsComm; }
void
ModelBase::setLocalNonCompositeWorldsComm(std::vector<WorldComm> const& _worldsComm) { M_localNonCompositeWorldsComm=_worldsComm; }
void
ModelBase::createWorldsComm() {}//warning

std::string const&
ModelBase::prefix() const { return M_prefix; }
std::string const&
ModelBase::subPrefix() const { return M_subPrefix; }

std::string const&
ModelBase::rootRepository() const { return M_modelRepository.root(); }

// verbose
bool
ModelBase::verbose() const { return M_verbose; }
bool
ModelBase::verboseAllProc() const { return M_verboseAllProc; }
void
ModelBase::log( std::string const& _className,std::string const& _functionName,std::string const& _msg ) const
{
    if (this->verbose()) FeelModels::Log( prefixvm(this->prefix(),_className),_functionName, _msg,
                                         this->worldComm(),this->verboseAllProc());
}

// info
std::string
ModelBase::filenameSaveInfo() const
{
    return M_filenameSaveInfo;
}
void
ModelBase::setFilenameSaveInfo( std::string const& s )
{
    M_filenameSaveInfo = s;
}
boost::shared_ptr<std::ostringstream>
ModelBase::getInfo() const
{
    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    return _ostr;
}
void
ModelBase::printInfo() const
{
    if ( this->verboseAllProc() )
        std::cout << this->getInfo()->str();
    else if (this->worldComm().isMasterRank() )
        std::cout << this->getInfo()->str();
}
void
ModelBase::saveInfo() const
{
    fs::path thepath = fs::path(this->rootRepository())/fs::path(this->filenameSaveInfo());
    if (this->worldComm().isMasterRank() )
    {
        std::ofstream file( thepath.string().c_str(), std::ios::out);
        file << this->getInfo()->str();
        file.close();
    }

    this->upload( thepath.string() );
}
void
ModelBase::printAndSaveInfo() const
{
    this->printInfo();
    this->saveInfo();
}

// timer
TimerToolBase &
ModelBase::timerTool(std::string const& key) const
{
    auto itFind = M_mapTimerTool.find( key );
    if ( itFind == M_mapTimerTool.end() )
    {
        this->addTimerTool(key,"applibase-defaultname.data");
        CHECK( M_mapTimerTool.find( key ) != M_mapTimerTool.end() ) << "key not exist";
        return *M_mapTimerTool[key];
    }
    else
        return *itFind->second;
}
void
ModelBase::addTimerTool(std::string const& key,std::string const& fileName) const
{
    CHECK( M_mapTimerTool.find( key ) == M_mapTimerTool.end() ) << "key already exist";
    if ( M_timersActivated )
    {
        auto myTimerTool = std::make_shared<TimerTool>(fileName,this->worldComm());
        myTimerTool->setReinitSaveFile( this->scalabilityReinitSaveFile() );
        myTimerTool->setSaveFileMasterRank( M_timersSaveFileMasterRank || M_timersSaveFileAll );
        myTimerTool->setSaveFileMax( M_timersSaveFileMax || M_timersSaveFileAll );
        myTimerTool->setSaveFileMin( M_timersSaveFileMin || M_timersSaveFileAll );
        myTimerTool->setSaveFileMean( M_timersSaveFileMean || M_timersSaveFileAll );
        M_mapTimerTool.emplace(std::make_pair( key,myTimerTool ) );
    }
    else
    {
        auto myTimerTool = std::make_shared<TimerToolNull>();
        M_mapTimerTool.emplace(std::make_pair( key,myTimerTool ) );
    }
}

// save assembly/solver scalability
bool
ModelBase::scalabilitySave() const { return M_scalabilitySave; }
bool
ModelBase::scalabilityReinitSaveFile() const { return M_scalabilityReinitSaveFile; }
void
ModelBase::setScalabilitySave( bool b )  { M_scalabilitySave=b; }
std::string
ModelBase::scalabilityPath() const { return M_scalabilityPath; }
void
ModelBase::setScalabilityPath( std::string const& s ) { M_scalabilityPath=s; }
std::string
ModelBase::scalabilityFilename() const { return M_scalabilityFilename; }
void
ModelBase::setScalabilityFilename( std::string const& s )  { M_scalabilityFilename=s; }

void
ModelBase::upload( std::string const& dataPath ) const
{
    if ( M_upload.isOperational() )
        M_upload.upload( dataPath );
}



} // namespace FeelModels

} // namespace Feel
