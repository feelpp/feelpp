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

namespace TabulateInformationTools
{
namespace FromJSON
{
tabulate_informations_ptr_t
tabulateInformationsModelFields( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    for ( auto const& el : jsonInfo.items() )
    {
        std::string const& fieldName = el.key();
        auto tabInfoField = TabulateInformationsSections::New( tabInfoProp );

        Feel::Table tabInfoBasic;
        TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoBasic, el.value(), tabInfoProp );
        tabInfoBasic.format()
            .setShowAllBorders( false )
            .setColumnSeparator(":")
            .setHasRowSeparator( false );
        tabInfoField->add( "", TabulateInformations::New( tabInfoBasic,tabInfoProp ) );
        if ( el.value().contains( "SymbolsExpr" ) )
        {
            auto const& jSE = el.value().at( "SymbolsExpr" );
            tabInfoField->add( "", TabulateInformationTools::FromJSON::tabulateInformationsSymbolsExpr( jSE,tabInfoProp.newByIncreasingVerboseLevel(),true) );
        }
        tabInfo->add( fieldName, tabInfoField );
    }
    return tabInfo;
}

} // namespace FromJSON
} // namespace TabulateInformationTools



namespace FeelModels {


void printToolboxApplication( std::string const& toolboxName, worldcomm_t const& worldComm )
{
    return;
    std::vector<std::string> all_lines;
    all_lines.push_back("███████╗███████╗███████╗██╗       ██╗         ██╗           ████████╗ ██████╗  ██████╗ ██╗     ██████╗  ██████╗ ██╗  ██╗███████╗███████╗");
    all_lines.push_back("██╔════╝██╔════╝██╔════╝██║       ██║         ██║           ╚══██╔══╝██╔═══██╗██╔═══██╗██║     ██╔══██╗██╔═══██╗╚██╗██╔╝██╔════╝██╔════╝");
    all_lines.push_back("█████╗  █████╗  █████╗  ██║   ██████████╗ ██████████╗          ██║   ██║   ██║██║   ██║██║     ██████╔╝██║   ██║ ╚███╔╝ █████╗  ███████╗");
    all_lines.push_back("██╔══╝  ██╔══╝  ██╔══╝  ██║   ╚═══██╔═══╝ ╚═══██╔═══╝  █████╗  ██║   ██║   ██║██║   ██║██║     ██╔══██╗██║   ██║ ██╔██╗ ██╔══╝  ╚════██║");
    all_lines.push_back("██║     ███████╗███████╗███████╗  ██║         ██║      ╚════╝  ██║   ╚██████╔╝╚██████╔╝███████╗██████╔╝╚██████╔╝██╔╝ ██╗███████╗███████║");
    all_lines.push_back("╚═╝     ╚══════╝╚══════╝╚══════╝  ╚═╝         ╚═╝              ╚═╝    ╚═════╝  ╚═════╝ ╚══════╝╚═════╝  ╚═════╝ ╚═╝  ╚═╝╚══════╝╚══════╝");

    std::vector<std::string> all_lines_app;
    if ( toolboxName == "cfpdes" )
    {
        all_lines_app.push_back(" ██████╗ ██████╗ ███████╗███████╗███████╗██╗ ██████╗██╗███████╗███╗   ██╗████████╗    ███████╗ ██████╗ ██████╗ ███╗   ███╗    ██████╗ ██████╗ ███████╗███████╗");
        all_lines_app.push_back("██╔════╝██╔═══██╗██╔════╝██╔════╝██╔════╝██║██╔════╝██║██╔════╝████╗  ██║╚══██╔══╝    ██╔════╝██╔═══██╗██╔══██╗████╗ ████║    ██╔══██╗██╔══██╗██╔════╝██╔════╝");
        all_lines_app.push_back("██║     ██║   ██║█████╗  █████╗  █████╗  ██║██║     ██║█████╗  ██╔██╗ ██║   ██║       █████╗  ██║   ██║██████╔╝██╔████╔██║    ██████╔╝██║  ██║█████╗  ███████╗");
        all_lines_app.push_back("██║     ██║   ██║██╔══╝  ██╔══╝  ██╔══╝  ██║██║     ██║██╔══╝  ██║╚██╗██║   ██║       ██╔══╝  ██║   ██║██╔══██╗██║╚██╔╝██║    ██╔═══╝ ██║  ██║██╔══╝  ╚════██║");
        all_lines_app.push_back("╚██████╗╚██████╔╝███████╗██║     ██║     ██║╚██████╗██║███████╗██║ ╚████║   ██║       ██║     ╚██████╔╝██║  ██║██║ ╚═╝ ██║    ██║     ██████╔╝███████╗███████║");
        all_lines_app.push_back(" ╚═════╝ ╚═════╝ ╚══════╝╚═╝     ╚═╝     ╚═╝ ╚═════╝╚═╝╚══════╝╚═╝  ╚═══╝   ╚═╝       ╚═╝      ╚═════╝ ╚═╝  ╚═╝╚═╝     ╚═╝    ╚═╝     ╚═════╝ ╚══════╝╚══════╝");
    }
    else if ( toolboxName == "heat" )
    {
        all_lines_app.push_back("██╗  ██╗███████╗ █████╗ ████████╗");
        all_lines_app.push_back("██║  ██║██╔════╝██╔══██╗╚══██╔══╝");
        all_lines_app.push_back("███████║█████╗  ███████║   ██║   ");
        all_lines_app.push_back("██╔══██║██╔══╝  ██╔══██║   ██║   ");
        all_lines_app.push_back("██║  ██║███████╗██║  ██║   ██║   ");
        all_lines_app.push_back("╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝   ╚═╝   ");
    }

    tabulate::Table tabInfo;
    // the a\b added is for fix the FontAlign::center, not very nice but work if we hide the left/right borders
    for ( int k=0;k<all_lines.size();++k )
        tabInfo.add_row({ "a\b" +all_lines[k] + "a\b" });
    for ( int k=0;k<all_lines_app.size();++k )
        tabInfo.add_row({ "a\b" +all_lines_app[k] +"a\b" });
    tabInfo.format()
        .hide_border()
        .multi_byte_characters(true)
        //.hide_border_top()
        .font_color(tabulate::Color::green/*red*/)
        .font_align(tabulate::FontAlign::center)
        //.font_background_color(tabulate::Color::green)
        ;

    if ( worldComm.isMasterRank() )
        std::cout << tabInfo << std::endl;
    worldComm.barrier();
}



namespace ToolboxesDetail
{
void removeTrailingSlash( std::string & s )
{
    if ( Feel::filename_is_dot( fs::path( s ).filename() ) )
        s = fs::path( s ).parent_path().string();
}
}

ModelBaseCommandLineOptions::ModelBaseCommandLineOptions( po::options_description const& _options )
{
    M_vm.emplace();
    auto mycmdparser = Environment::commandLineParser();
    po::parsed_options parsed = mycmdparser.options( _options ).
        style(po::command_line_style::allow_long | po::command_line_style::long_allow_adjacent | po::command_line_style::long_allow_next).
        allow_unregistered().run();
    po::store(parsed,*M_vm);
    for ( auto & configFile : Environment::configFiles() )
    {
        std::istringstream & iss = std::get<1>( configFile );
        po::store(po::parse_config_file(iss, _options,true), *M_vm);
    }
    po::notify(*M_vm);
}


ModelBaseRepository::ModelBaseRepository( std::string const& rootDirWithoutNumProc, bool use_npSubDir, std::string const& exprRepository )
{
    if ( rootDirWithoutNumProc.empty() )
    {
        M_rootRepositoryWithoutNumProc = Environment::appRepositoryWithoutNumProc();
        if ( use_npSubDir )
            M_rootRepositoryWithNumProc = Environment::appRepository();
        else
            M_rootRepositoryWithNumProc = M_rootRepositoryWithoutNumProc;
        M_exprRepository = Environment::exprRepository();
    }
    else
    {
        M_rootRepositoryWithoutNumProc = Environment::expand( rootDirWithoutNumProc );
        if ( fs::path( M_rootRepositoryWithoutNumProc ).is_relative() )
            M_rootRepositoryWithoutNumProc = (fs::path(Environment::rootRepository())/fs::path(M_rootRepositoryWithoutNumProc)).string();

        if ( use_npSubDir )
        {
            std::string npSubDir = (boost::format( "np_%1%" ) % Environment::worldComm().localSize() ).str();
            M_rootRepositoryWithNumProc = ( fs::path(M_rootRepositoryWithoutNumProc) / fs::path( npSubDir ) ).string();
        }
        else
            M_rootRepositoryWithNumProc = M_rootRepositoryWithoutNumProc;
        M_exprRepository = ( fs::path(M_rootRepositoryWithoutNumProc) / fs::path( "exprs" ) ).string();
    }

    if ( !exprRepository.empty() )
    {
        M_exprRepository = Environment::expand( exprRepository );
        if ( fs::path( M_exprRepository ).is_relative() )
            M_exprRepository = (fs::path(M_rootRepositoryWithoutNumProc)/fs::path( M_exprRepository )).string();
    }

    ToolboxesDetail::removeTrailingSlash( M_rootRepositoryWithoutNumProc );
    ToolboxesDetail::removeTrailingSlash( M_rootRepositoryWithNumProc );
    ToolboxesDetail::removeTrailingSlash( M_exprRepository );
}

ModelBaseUpload::ModelBaseUpload( std::string const& desc, std::string const& basePath, worldcomm_ptr_t const& worldComm )
    :
    M_basePath( basePath )
{
    if ( desc.empty() )
        return;
    M_remoteData.reset( new RemoteData( desc, worldComm ) );
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


ModelBase::ModelBase( std::string const& prefix, std::string const& keyword,
                      worldcomm_ptr_t const& worldComm,
                      std::string const& subPrefix,
                      ModelBaseRepository const& modelRep,
                      ModelBaseCommandLineOptions const& modelCmdLineOpt )
    :
    super_type( "Toolboxes",keyword/*, prefix*/ ),
    M_worldComm(worldComm),
    M_worldsComm( {worldComm} ),
    M_localNonCompositeWorldsComm( { worldComm } ),
    M_modelCommandLineOptions( modelCmdLineOpt ),
    M_prefix( prefix ),
    M_subPrefix( subPrefix ),
    M_keyword( keyword ),
    M_modelRepository( modelRep ),
    M_verbose( boption(_name="verbose",_prefix=this->prefix(),_vm=this->clovm()) ),
    M_verboseAllProc( boption(_name="verbose_allproc",_prefix=this->prefix(),_vm=this->clovm()) ),
    M_timersActivated( boption(_name="timers.activated",_prefix=this->prefix(),_vm=this->clovm()) ),
    M_timersSaveFileMasterRank( boption(_name="timers.save-master-rank",_prefix=this->prefix(),_vm=this->clovm()) ),
    M_timersSaveFileMax( boption(_name="timers.save-max",_prefix=this->prefix(),_vm=this->clovm()) ),
    M_timersSaveFileMin( boption(_name="timers.save-min",_prefix=this->prefix(),_vm=this->clovm()) ),
    M_timersSaveFileMean( boption(_name="timers.save-mean",_prefix=this->prefix(),_vm=this->clovm()) ),
    M_timersSaveFileAll( boption(_name="timers.save-all",_prefix=this->prefix(),_vm=this->clovm()) ),
    M_scalabilitySave( boption(_name="scalability-save",_prefix=this->prefix(),_vm=this->clovm()) ),
    M_scalabilityReinitSaveFile( boption(_name="scalability-reinit-savefile",_prefix=this->prefix(),_vm=this->clovm()) ),
    M_isUpdatedForUse( false ),
    M_upload( soption(_name="upload",_prefix=this->prefix(),_vm=this->clovm()), this->repository().rootWithoutNumProc(), M_worldComm )
{
    if (this->clovm().count(prefixvm(this->prefix(),"scalability-path")))
        M_scalabilityPath = soption(_name="scalability-path",_prefix=this->prefix(),_vm=this->clovm());
    else
        M_scalabilityPath = this->repository().rootWithoutNumProc();

    if (this->clovm().count(prefixvm(this->prefix(),"scalability-filename")))
        M_scalabilityFilename = soption(_name="scalability-filename",_prefix=this->prefix(),_vm=this->clovm());
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

ModelBase::worldcomm_ptr_t  const&
ModelBase::worldCommPtr() const { return M_worldComm; }
ModelBase::worldcomm_ptr_t  &
ModelBase::worldCommPtr() { return M_worldComm; }
WorldComm &
ModelBase::worldComm() { return *M_worldComm; }
WorldComm const&
ModelBase::worldComm() const { return *M_worldComm; }
ModelBase::worldscomm_ptr_t &
ModelBase::worldsComm()  { return M_worldsComm; }
ModelBase::worldscomm_ptr_t const&
ModelBase::worldsComm() const { return M_worldsComm; }
void
ModelBase::setWorldsComm( worldscomm_ptr_t & _worldsComm) { M_worldsComm=_worldsComm; }
ModelBase::worldscomm_ptr_t &
ModelBase::localNonCompositeWorldsComm() { return M_localNonCompositeWorldsComm; }
ModelBase::worldscomm_ptr_t const&
ModelBase::localNonCompositeWorldsComm() const { return M_localNonCompositeWorldsComm; }
void
ModelBase::setLocalNonCompositeWorldsComm( worldscomm_ptr_t& _worldsComm) { M_localNonCompositeWorldsComm=_worldsComm; }
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

void
ModelBase::updateInformationObject( nl::json & p ) const
{
    if ( p.contains( "prefix" ) )
        return;
    p.emplace( "prefix", this->prefix() );
    p.emplace( "keyword", this->keyword() );
    p.emplace( "root repository", this->rootRepository() );
    p.emplace( "expr repository", this->repository().expr() );
    p.emplace( "number of processus", this->worldComm().localSize() );
}

tabulate_informations_ptr_t
ModelBase::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    Feel::Table tabInfo;
    TabulateInformationTools::FromJSON::addKeyToValues( tabInfo, jsonInfo, tabInfoProp, { "prefix","keyword","root repository","expr eepository", "number of processus" } );
    //TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfo, jsonInfo, tabInfoProp ); // bad ordering due to boost properties
    tabInfo.format()
        .setShowAllBorders( false )
        .setColumnSeparator(":")
        .setHasRowSeparator( false );
    return TabulateInformations::New( tabInfo, tabInfoProp );
}

tabulate_informations_ptr_t
ModelBase::tabulateInformations() const
{
    Environment::journalCheckpoint();
#if 0
    pt::ptree pt;
    this->updateInformationObject( pt );
    std::ostringstream pt_ostr;
    write_json( pt_ostr, pt );
    std::istringstream pt_istream( pt_ostr.str() );
    nl::json jsonInfo;
    pt_istream >> jsonInfo;
#else
    nl::json jsonInfo;
    this->updateInformationObject( jsonInfo );
#endif

    auto tabInfo = this->tabulateInformations( jsonInfo, TabulateInformationProperties{} );

    auto tabRes = TabulateInformationsSections::New();
    std::string title = (boost::format("Toolbox : %1%")%this->keyword()).str();
    tabRes->add( title, tabInfo );

    return tabRes;
}

std::shared_ptr<std::ostringstream>
ModelBase::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    return _ostr;
}
void
ModelBase::printInfo( tabulate_informations_ptr_t const& tabInfos ) const
{
    if ( this->worldComm().isMasterRank() )
    {
        std::cout << this->getInfo()->str();
        std::cout << *tabInfos << std::endl;
    }
}
void
ModelBase::saveInfo( tabulate_informations_ptr_t const& tabInfos ) const
{
    std::string filename_ascii = prefixvm(this->keyword(),"informations.txt");
    std::string filename_adoc = prefixvm(this->keyword(),"informations.adoc");
    std::string filepath_ascii = (fs::path(this->rootRepository())/filename_ascii).string();
    std::string filepath_adoc = (fs::path(this->rootRepository())/filename_adoc).string();
    if ( this->worldComm().isMasterRank() )
    {
        std::ofstream file_ascii( filepath_ascii, std::ios::out);
        file_ascii << this->getInfo()->str();
        file_ascii.close();

        std::ofstream file_adoc( filepath_adoc, std::ios::out);
        file_adoc << ":sectnums:" << "\n";
        file_adoc << tabInfos->exporterAsciiDoc() << "\n";
        file_adoc.close();
    }

    this->upload( filepath_ascii );
    this->upload( filepath_adoc );
}
void
ModelBase::printAndSaveInfo() const
{
    auto tabInfo = this->tabulateInformations();
    this->printInfo( tabInfo );
    this->saveInfo( tabInfo );
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
        auto myTimerTool = std::make_shared<TimerTool>(fileName,this->worldCommPtr());
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
