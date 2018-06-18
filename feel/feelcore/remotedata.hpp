/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannnes@feelpp.org>
 Date: 2 April. 2018

 Copyright (C) 2018 Feel++ Consortium

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#ifndef FEELPP_CORE_REMOTEDATA_HPP
#define FEELPP_CORE_REMOTEDATA_HPP 1

#include <boost/property_tree/json_parser.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>

namespace Feel
{

/**
 * \brief class which manage downloads from an url or github and download/upload with the Girder database by using the libcurl
 */
struct RemoteData
{
    RemoteData( std::string const& desc, WorldComm const& worldComm = Environment::worldComm() );
    RemoteData( RemoteData const& ) = default;
    RemoteData( RemoteData && ) = default;

    //! return true if enough information are available for download a file
    bool canDownload() const;

    //! return true if enough information are available for upload data
    bool canUpload() const;

    //! download the file
    //! @param dir : the directory where the file is downloaded
    //! @param filename : the filename of the downloaded file
    //! @return : the path of the downloaded file
    std::vector<std::string> download( std::string const& dir = Environment::downloadsRepository(), std::string const& filename = "" ) const;

    //! upload data on a remote storage
    //! @param dataPath : a path of a file or a folder
    //! @param parentId : id where folder are created, empty say to use folder id in the desc
    void upload( std::string const& dataPath, std::string const& parentId = "" ) const;

    //! create folders hierarchy on remote storage
    //! @param folderPath : the folder hierarchy
    //! @param parentId : id where folder are created, empty say to use folder id in the desc
    //! @return : vector of (subdir name, subdir id)
    std::vector<std::pair<std::string,std::string>>
    createFolder( std::string const& folderPath, std::string const& parentId = "" ) const;

    class URL
    {
    public :
        URL( std::string const& url, WorldComm const& worldComm = Environment::worldComm() );
        URL( URL const& ) = default;
        URL( URL && ) = default;

        //! return true if the url is valid
        bool isValid() const;

        //! download a file from the url
        //! @param dir : the directory where the file is downloaded
        //! @param filename : the filename of the downloaded file
        //! @return : the path of the downloaded file
        std::string download( std::string const& dir = Environment::downloadsRepository(), std::string const& filename = "" ) const;
    private :
        WorldComm const& M_worldComm;
        std::string M_url;
        std::string M_protocol, M_domain, M_port, M_path, M_query;
    };

    class Github
    {
    public :
        //! init from a description :
        //!   github:{owner:feelpp,repo:feelpp,branch:develop,path:toolboxes/fluid/TurekHron,token:xxxxx}
        Github( std::string const& desc, WorldComm const& worldComm = Environment::worldComm() );
        Github( Github const& ) = default;
        Github( Github && ) = default;

        //! return true if the github is initialized from a desc
        bool isInit() const;

        //! download file/folder from the github desc
        //! @param dir : the directory where the file is downloaded
        //! @return : vector of path of the downloaded file or the path of downloaded folder
        std::vector<std::string> download( std::string const& dir = Environment::downloadsRepository() ) const;
    private :
        void downloadFolderRecursively( pt::ptree const& ptree, std::string const& dir ) const;

    private :
        WorldComm const& M_worldComm;
        std::string M_owner, M_repo, M_branch, M_path, M_token;
    };

    class Girder
    {
    public :
        //! init from a description : github:{owner:feelpp,repo:feelpp,branch:develop,path:toolboxes/fluid/TurekHron,token:xxxxx}
        //!   girder:{url:https://girder.math.unistra.fr,file:5ac722e9b0e9574027047886,token:xxxxx}
        //!   girder:{url:https://girder.math.unistra.fr,file:[5ac7253ab0e957402704788d,5ac722e9b0e9574027047886],token:xxxxx}
        Girder( std::string const& desc, WorldComm const& worldComm = Environment::worldComm() );
        Girder( Girder const& ) = default;
        Girder( Girder && ) = default;

        //! set folder ids to only one folder id
        void setFolderIds( std::string const& folderId );

        //! return true if the github is initialized from a desc
        bool isInit() const;

        //! return true if enough information are available for download a file
        bool canDownload() const;

        //! return true if enough information are available for upload data
        bool canUpload() const;

        //! download file/folder from the girder desc
        //! @param dir : the directory where the file is downloaded
        //! @return : vector of path of the downloaded file or the path of downloaded folder
        std::vector<std::string> download( std::string const& dir = Environment::downloadsRepository() ) const;

        //! upload data on Girder
        //! @param dataPath : a path of a file or a folder
        //! @param parentId : id where folder are created, empty say to use folder id in the desc
        void upload( std::string const& dataPath, std::string const& parentId = "" ) const;

        //! create folders hierarchy on Girder
        //! @param folderPath : the folder hierarchy
        //! @param parentId : id where folder are created, empty say to use folder id in the desc
        //! @return : vector of (subdir name, subdir id)
        std::vector<std::pair<std::string,std::string>>
        createFolder( std::string const& folderPath, std::string const& parentId = "" ) const;

    private :
        std::string downloadFile( std::string const& fileId, std::string const& dir ) const;
        void uploadRecursively( std::string const& dataPath, std::string const& parentId, std::string const& token ) const;
        void uploadFile( std::string const& filePath, std::string const& parentId, std::string const& token ) const;
        std::string createFolderImpl( std::string const& folderName, std::string const& parentId, std::string const& token ) const;
        std::string createToken( int duration = 1 ) const;
        void removeToken( std::string const& token ) const;

    private :
        WorldComm const& M_worldComm;
        std::string M_url, M_apiKey, M_token;
        std::set<std::string> M_fileIds, M_folderIds;//, M_itemId;
    };

private :
    boost::optional<URL> M_url;
    boost::optional<Github> M_github;
    boost::optional<Girder> M_girder;
};

}

#endif
