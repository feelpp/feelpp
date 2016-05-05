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

namespace Feel
{

namespace FeelModels
{

ModelBase::ModelBase( std::string const& prefix,
                      WorldComm const& worldComm,
                      std::string const& subPrefix,
                      std::string const& rootRepository )
    : M_worldComm( worldComm ),
      M_worldsComm( std::vector<WorldComm>( 1, worldComm ) ),
      M_localNonCompositeWorldsComm( std::vector<WorldComm>( 1, worldComm ) ),
      M_prefix( prefix ),
      M_subPrefix( subPrefix ),
      M_rootRepositoryWithoutNumProc( rootRepository ),
      M_verbose( boption( _name = "verbose", _prefix = this->prefix() ) ),
      M_verboseAllProc( boption( _name = "verbose_allproc", _prefix = this->prefix() ) ),
      M_filenameSaveInfo( prefixvm( this->prefix(), prefixvm( this->subPrefix(), "appli.info" ) ) ),
      M_timersActivated( boption( _name = "timers.activated", _prefix = this->prefix() ) ),
      M_timersSaveFileMasterRank( boption( _name = "timers.save-master-rank", _prefix = this->prefix() ) ),
      M_timersSaveFileMax( boption( _name = "timers.save-max", _prefix = this->prefix() ) ),
      M_timersSaveFileMin( boption( _name = "timers.save-min", _prefix = this->prefix() ) ),
      M_timersSaveFileMean( boption( _name = "timers.save-mean", _prefix = this->prefix() ) ),
      M_timersSaveFileAll( boption( _name = "timers.save-all", _prefix = this->prefix() ) ),
      M_scalabilitySave( boption( _name = "scalability-save", _prefix = this->prefix() ) ),
      M_scalabilityReinitSaveFile( boption( _name = "scalability-reinit-savefile", _prefix = this->prefix() ) )
{
    if ( M_rootRepositoryWithoutNumProc.empty() )
        M_rootRepositoryWithoutNumProc = ModelBase::rootRepositoryByDefault();
    M_rootRepositoryWithoutNumProc = Environment::expand( M_rootRepositoryWithoutNumProc );
    if ( fs::path( M_rootRepositoryWithoutNumProc ).is_relative() )
        M_rootRepositoryWithoutNumProc = ( fs::path( Environment::rootRepository() ) / fs::path( M_rootRepositoryWithoutNumProc ) ).string();
    std::string npSubDir = ( boost::format( "np_%1%" ) % this->worldComm().localSize() ).str();
    M_rootRepositoryWithNumProc = ( fs::path( M_rootRepositoryWithoutNumProc ) / fs::path( npSubDir ) ).string();

    if ( Environment::vm().count( prefixvm( this->prefix(), "scalability-path" ) ) )
        M_scalabilityPath = Environment::vm()[prefixvm( this->prefix(), "scalability-path" )].as<std::string>();
    else
        M_scalabilityPath = this->rootRepositoryWithoutNumProc();

    if ( Environment::vm().count( prefixvm( this->prefix(), "scalability-filename" ) ) )
        M_scalabilityFilename = Environment::vm()[prefixvm( this->prefix(), "scalability-filename" )].as<std::string>();
    else
        M_scalabilityFilename = this->prefix() + ".scalibility";
}
ModelBase::~ModelBase()
{
}

WorldComm const&
ModelBase::worldComm() const
{
    return M_worldComm;
}
std::vector<WorldComm> const&
ModelBase::worldsComm() const
{
    return M_worldsComm;
}
void ModelBase::setWorldsComm( std::vector<WorldComm> const& _worldsComm )
{
    M_worldsComm = _worldsComm;
}
std::vector<WorldComm> const&
ModelBase::localNonCompositeWorldsComm() const
{
    return M_localNonCompositeWorldsComm;
}
void ModelBase::setLocalNonCompositeWorldsComm( std::vector<WorldComm> const& _worldsComm )
{
    M_localNonCompositeWorldsComm = _worldsComm;
}
void ModelBase::createWorldsComm()
{
} //warning

std::string const&
ModelBase::prefix() const
{
    return M_prefix;
}
std::string const&
ModelBase::subPrefix() const
{
    return M_subPrefix;
}

po::variables_map const&
ModelBase::vm() const
{
    return Environment::vm();
}

std::string const&
ModelBase::rootRepository() const
{
    return this->rootRepositoryWithNumProc();
}
std::string const&
ModelBase::rootRepositoryWithoutNumProc() const
{
    return M_rootRepositoryWithoutNumProc;
}
std::string const&
ModelBase::rootRepositoryWithNumProc() const
{
    return M_rootRepositoryWithNumProc;
}

std::string
ModelBase::rootRepositoryByDefault()
{
    return soption( _name = "exporter.directory" );
}

// verbose
bool ModelBase::verbose() const
{
    return M_verbose;
}
bool ModelBase::verboseAllProc() const
{
    return M_verboseAllProc;
}
void ModelBase::log( std::string const& _className, std::string const& _functionName, std::string const& _msg ) const
{
    if ( this->verbose() ) FeelModels::Log( prefixvm( this->prefix(), _className ), _functionName, _msg,
                                            this->worldComm(), this->verboseAllProc() );
}

// info
std::string
ModelBase::filenameSaveInfo() const
{
    return M_filenameSaveInfo;
}
void ModelBase::setFilenameSaveInfo( std::string const& s )
{
    M_filenameSaveInfo = s;
}
boost::shared_ptr<std::ostringstream>
ModelBase::getInfo() const
{
    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    return _ostr;
}
void ModelBase::printInfo() const
{
    if ( this->verboseAllProc() )
        std::cout << this->getInfo()->str();
    else if ( this->worldComm().isMasterRank() )
        std::cout << this->getInfo()->str();
}
void ModelBase::saveInfo() const
{
    if ( this->worldComm().isMasterRank() )
    {
        fs::path thepath = fs::path( this->rootRepository() ) / fs::path( this->filenameSaveInfo() );
        std::ofstream file( thepath.string().c_str(), std::ios::out );
        file << this->getInfo()->str();
        file.close();
    }
}
void ModelBase::printAndSaveInfo() const
{
    this->printInfo();
    this->saveInfo();
}

// timer
TimerToolBase&
ModelBase::timerTool( std::string const& key ) const
{
    auto itFind = M_mapTimerTool.find( key );
    if ( itFind == M_mapTimerTool.end() )
    {
        this->addTimerTool( key, "applibase-defaultname.data" );
        CHECK( M_mapTimerTool.find( key ) != M_mapTimerTool.end() ) << "key not exist";
        return *M_mapTimerTool[key];
    }
    else
        return *itFind->second;
}
void ModelBase::addTimerTool( std::string const& key, std::string const& fileName ) const
{
    CHECK( M_mapTimerTool.find( key ) == M_mapTimerTool.end() ) << "key already exist";
    if ( M_timersActivated )
    {
        auto myTimerTool = std::make_shared<TimerTool>( fileName, this->worldComm() );
        myTimerTool->setReinitSaveFile( this->scalabilityReinitSaveFile() );
        myTimerTool->setSaveFileMasterRank( M_timersSaveFileMasterRank || M_timersSaveFileAll );
        myTimerTool->setSaveFileMax( M_timersSaveFileMax || M_timersSaveFileAll );
        myTimerTool->setSaveFileMin( M_timersSaveFileMin || M_timersSaveFileAll );
        myTimerTool->setSaveFileMean( M_timersSaveFileMean || M_timersSaveFileAll );
        M_mapTimerTool.emplace( std::make_pair( key, myTimerTool ) );
    }
    else
    {
        auto myTimerTool = std::make_shared<TimerToolNull>();
        M_mapTimerTool.emplace( std::make_pair( key, myTimerTool ) );
    }
}

// save assembly/solver scalability
bool ModelBase::scalabilitySave() const
{
    return M_scalabilitySave;
}
bool ModelBase::scalabilityReinitSaveFile() const
{
    return M_scalabilityReinitSaveFile;
}
void ModelBase::setScalabilitySave( bool b )
{
    M_scalabilitySave = b;
}
std::string
ModelBase::scalabilityPath() const
{
    return M_scalabilityPath;
}
void ModelBase::setScalabilityPath( std::string const& s )
{
    M_scalabilityPath = s;
}
std::string
ModelBase::scalabilityFilename() const
{
    return M_scalabilityFilename;
}
void ModelBase::setScalabilityFilename( std::string const& s )
{
    M_scalabilityFilename = s;
}

} // namespace FeelModels

} // namespace Feel
