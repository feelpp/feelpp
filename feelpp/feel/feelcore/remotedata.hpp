/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannnes@feelpp.org>
 Date: 2 April. 2018

 Copyright (C) 2018-2024 Feel++ Consortium

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

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/json.hpp>

namespace Feel
{
using namespace Feel;

/**
 * \brief Class which manages downloads from a URL or GitHub and download/upload with the Girder database by using libcurl
 */
struct RemoteData
{
    struct FolderInfo;
    struct ItemInfo;
    struct FileInfo;

    RemoteData( std::string const& desc, worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );
    FEELPP_DEPRECATED RemoteData( std::string const& desc, WorldComm const& worldComm )
        : RemoteData( desc, const_cast<WorldComm&>( worldComm ).shared_from_this() ) {}
    RemoteData( RemoteData const& ) = default;
    RemoteData( RemoteData&& ) = default;
    ~RemoteData() = default;

    //! Return world comm
    WorldComm& worldComm() const { return *M_worldComm; }

    //! Set the WorldComm
    void setWorldComm( WorldComm& wc ) { M_worldComm = wc.shared_from_this(); }

    //! Return true if enough information is available to download a file
    bool canDownload() const;

    //! Return true if enough information is available to upload data
    bool canUpload() const;

    //! Download the file
    //! @param dir : the directory where the file is downloaded
    //! @param filename : the filename of the downloaded file
    //! @return : the path of the downloaded file
    std::vector<std::string> download( std::string const& dir = Environment::downloadsRepository(), std::string const& filename = "" ) const;

    //! Upload data on a remote storage
    //! @param dataPath : a path of a file or a folder
    //! @param parentId : id where folder is created, empty means to use folder id in the desc
    //! @param sync : apply MPI synchronization with returned infos (else only master rank has these infos)
    //! @return : vector of paths of the uploaded files or folders
    std::vector<std::string>
    upload( std::string const& dataPath, std::string const& parentId = "", bool sync = true ) const;

    //! Upload data on Girder
    //! @param dataToUpload : vector of (paths of a file or a folder, ids where folder is created, empty means to use folder id in the desc)
    //! @param sync : apply MPI synchronization with returned infos (else only master rank has these infos)
    //! @return : vector of vectors of file ids uploaded
    std::vector<std::vector<std::string>>
    upload( std::vector<std::pair<std::string, std::string>> const& dataToUpload, bool sync = true ) const;

    //! Replace contents of a file
    //! @param filePath : path of new file
    //! @param fileId : id of the file to replace
    void replaceFile( std::string const& filePath, std::string const& fileId ) const;

    //! Replace contents of files
    //! @param filesToReplace : vector of (path of new file, file id to replace)
    void replaceFile( std::vector<std::pair<std::string, std::string>> const& filesToReplace ) const;

    //! Create folders hierarchy on remote storage
    //! @param folderPath : the folder hierarchy
    //! @param parentId : id where folder is created, empty means to use folder id in the desc
    //! @param sync : apply MPI synchronization with returned infos (else only master rank has these infos)
    //! @return : vector of (subdir name, subdir id)
    std::vector<std::pair<std::string, std::string>>
    createFolder( std::string const& folderPath, std::string const& parentId = "", bool sync = true ) const;

    //! Create item  on remote storage
    //! @param itemPath : the item 
    //! @param parentId : id where item is created, empty means to use folder id in the desc
    //! @param sync : apply MPI synchronization with returned infos (else only master rank has these infos)
    std::vector<std::pair<std::string, std::string>>
    createItem( std::string const& itemPath, std::string const& parentId, bool sync = true ) const;

    //! Content info data structure
    class ContentsInfo : public std::tuple<std::vector<std::shared_ptr<FolderInfo>>, std::vector<std::shared_ptr<ItemInfo>>, std::vector<std::shared_ptr<FileInfo>>>
    {
      public:
        using base = std::tuple<std::vector<std::shared_ptr<FolderInfo>>, std::vector<std::shared_ptr<ItemInfo>>, std::vector<std::shared_ptr<FileInfo>>>;
        ContentsInfo() = default;
        explicit ContentsInfo( base const& b )
            : base( b ) {}
        ~ContentsInfo() = default;

        std::vector<std::shared_ptr<FolderInfo>> folderInfo() { return std::get<0>( *this ); }
        std::vector<std::shared_ptr<ItemInfo>> itemInfo() { return std::get<1>( *this ); }
        std::vector<std::shared_ptr<FileInfo>> fileInfo() { return std::get<2>( *this ); }
    };

    //! Lookup a resource by path
    //! @param path : the path to the resource 
    //! @param token : authentication token
    //! @return : JSON object containing resource information
    nl::json resourceLookup( const std::string& path, const std::string& token = "" ) const;

    //! Delete a resource by id
    //! @param resourceId : the id of the resource in Girder
    //! @param token : authentication token
    //! @return : true if the resource was deleted
    bool deleteResource( const nl::json& resourceId, const std::string& token = "" ) const;

    //! Get contents of remote data (folder, item, file)
    //! @return : (Folders info, Items info, Files info)
    ContentsInfo contents() const;

    class URL
    {
      public:
        URL( std::string const& url, WorldComm& worldComm = Environment::worldComm() );
        URL( URL const& ) = default;
        URL( URL&& ) = default;

        //! Return true if the URL is valid
        bool isValid() const;

        //! Download a file from the URL
        //! @param dir : the directory where the file is downloaded
        //! @param filename : the filename of the downloaded file
        //! @return : the path of the downloaded file
        std::string download( std::string const& dir = Environment::downloadsRepository(), std::string const& filename = "" ) const;

      private:
        std::shared_ptr<WorldComm> M_worldComm;
        std::string M_url;
        std::string M_protocol, M_domain, M_port, M_path, M_query;
    };

    class Github
    {
      public:
        //! Initialize from a description:
        //!   github:{owner:feelpp,repo:feelpp,branch:develop,path:toolboxes/fluid/TurekHron,token:xxxxx}
        Github( std::string const& desc, WorldComm& worldComm = Environment::worldComm() );
        Github( Github const& ) = default;
        Github( Github&& ) = default;

        //! Return true if the GitHub is initialized from a desc
        bool isInit() const;

        //! Download file/folder from the GitHub desc
        //! @param dir : the directory where the file is downloaded
        //! @return : vector of paths of the downloaded files or the path of downloaded folder
        std::vector<std::string> download( std::string const& dir = Environment::downloadsRepository() ) const;

      private:
        std::vector<std::string> downloadImpl( std::string const& dir ) const;
        std::tuple<bool, std::string> downloadFolderRecursively( nl::json const& jsonResponse, std::string const& dir ) const;

        static std::string errorMessage( nl::json const& jsonResponse, std::string const& defaultMsg = "", uint16_type statusCode = invalid_uint16_type_value );

      private:
        std::shared_ptr<WorldComm> M_worldComm;
        std::string M_owner, M_repo, M_branch, M_path, M_token;
    };

    class Girder
    {
      public:
        //! Initialize from a description:
        //!   girder:{url:https://girder.math.unistra.fr,file:5ac722e9b0e9574027047886,token:xxxxx}
        //!   girder:{url:https://girder.math.unistra.fr,file:[5ac7253ab0e957402704788d,5ac722e9b0e9574027047886],token:xxxxx}
        Girder( std::string const& desc, WorldComm& worldComm = Environment::worldComm() );
        Girder( Girder const& ) = default;
        Girder( Girder&& ) = default;

        //! Set folder ids to only one folder id
        void setFolderIds( std::string const& folderId );

        //! Return true if the Girder remote server is defined from a desc
        bool isInit() const;

        //! Return true if enough information is available to download a file
        bool canDownload() const;

        //! Return true if enough information is available to upload data
        bool canUpload() const;

        //! Download file/folder from the Girder desc
        //! @param dir : the directory where the file is downloaded
        //! @return : vector of paths of the downloaded files or the path of downloaded folder
        std::vector<std::string> download( std::string const& dir = Environment::downloadsRepository() ) const;

        //! Download file/folder/item from the Girder desc
        //! @param dir : the directory where the file is downloaded
        //! @param path : the path of the file/folder/item
        //! @return : vector of paths of the downloaded files or the path of downloaded folder
        std::vector<std::string> download( const std::string& dir, const std::string& path ) const;

        //! Upload data on Girder
        //! @param dataPath : a path of a file or a folder
        //! @param parentId : id where folder is created, empty means to use folder id in the desc
        //! @param sync : apply MPI synchronization with returned infos (else only master rank has these infos)
        //! @return : vector of file ids uploaded
        std::vector<std::string>
        upload( std::string const& dataPath, std::string const& parentId = "", bool sync = true ) const;

        //! Upload data on Girder
        //! @param dataToUpload : vector of (paths of a file or a folder, ids where folder is created, empty means to use folder id in the desc)
        //! @param sync : apply MPI synchronization with returned infos (else only master rank has these infos)
        //! @return : vector of vectors of file ids uploaded
        std::vector<std::vector<std::string>>
        upload( std::vector<std::pair<std::string, std::string>> const& dataToUpload, bool sync = true ) const;


        //! Replace contents of a file
        //! @param filePath : path of new file
        //! @param fileId : id of the file to replace
        void replaceFile( std::string const& filePath, std::string const& fileId ) const;

        //! Replace contents of files
        //! @param filesToReplace : vector of (path of new file, file id to replace)
        void replaceFile( std::vector<std::pair<std::string, std::string>> const& filesToReplace ) const;

        //! Create folders hierarchy on Girder
        //! @param folderPath : the folder hierarchy
        //! @param parentId : id where folder is created, empty means to use folder id in the desc
        //! @param sync : apply MPI synchronization with returned infos (else only master rank has these infos)
        //! @return : vector of (subdir name, subdir id)
        std::vector<std::pair<std::string, std::string>>
        createFolder( std::string const& folderPath, std::string const& parentId = "", bool sync = true ) const;

        //! Create item  on remote storage
        //! @param itemPath : the item 
        //! @param parentId : id where item is created, empty means to use folder id in the desc
        //! @param sync : apply MPI synchronization with returned infos (else only master rank has these infos)
        std::vector<std::pair<std::string, std::string>>
        createItem( std::string const& itemPath, std::string const& parentId, bool sync = true ) const;

        //! Get contents of remote data (folder, item, file)
        //! @return : (Folders info, Items info, Files info)
        std::tuple<std::vector<std::shared_ptr<FolderInfo>>, std::vector<std::shared_ptr<ItemInfo>>, std::vector<std::shared_ptr<FileInfo>>>
        contents() const;

        //! Lookup a resource by path
        //! @param path : the path to the resource in Girder
        //! @param token : authentication token
        //! @return : JSON object containing resource information
        nl::json resourceLookup( const std::string& path, const std::string& token = "" ) const;

        //! Delete a resource by id
        //! @param resourceId : the id of the resource in Girder
        //! @param token : authentication token
        //! @return : true if the resource was deleted
        bool deleteResource( const nl::json& resourceId, const std::string& token = "" ) const;

      private:
        std::string downloadFile( std::string const& fileId, std::string const& dir, std::string const& token ) const;
        std::string downloadFolder( std::string const& folderId, std::string const& dir, std::string const& token ) const;
        std::vector<std::string> uploadRecursively( std::string const& dataPath, std::string const& parentId, std::string const& token ) const;
        //std::string uploadFileImpl( std::string const& filePath, std::string const& parentId, std::string const& token ) const;
        std::string uploadFileImpl(const std::string& filepath, const std::string& parentId, const std::string& token, const std::string& parentType) const;

        void uploadDirectoryToFolder(const std::string& localDir, const std::string& parentFolderId, const std::string& token, std::vector<std::string>& uploadedResources) const;
        void uploadFilesToItem(const std::string& localPath, const std::string& itemId, const std::string& token, std::vector<std::string>& uploadedResources) const;
        void replaceFileImpl( std::string const& filePath, std::string const& fileId, std::string const& token ) const;
        //std::string createFolderImpl( std::string const& folderName, std::string const& parentId, std::string const& token ) const;
        std::string createToken( int duration = 1 ) const;
        void removeToken( std::string const& token ) const;
        std::string createItemImpl(const std::string& itemName, const std::string& parentFolderId, const std::string& token) const;
        std::string createFolderImpl(const std::string& folderName, const std::string& parentId, const std::string& token, const std::string& parentType = "folder") const;
        std::string initializeUpload(const std::string& filename, std::streamsize fileSize, const std::string& parentId, const std::string& parentType, const std::string& token) const;
        std::string uploadChunk(const std::string& uploadId, std::streamsize offset, const char* data, std::streamsize size, const std::string& token, bool isFinalChunk = false) const;


        std::shared_ptr<FileInfo>
        fileInfoImpl( std::string const& fileId, std::string const& token ) const;
        std::shared_ptr<FolderInfo>
        folderInfoImpl( std::string const& folderId, std::string const& token ) const;
        std::shared_ptr<ItemInfo>
        itemInfoImpl( std::string const& itemId, std::string const& token ) const;
        std::shared_ptr<FolderInfo>
        folderContentsImpl( std::string const& folderId, std::string const& token ) const;
        void updateFilesImpl( std::shared_ptr<RemoteData::ItemInfo> itemInfo, std::string const& token ) const;

        static std::string errorMessage( nl::json const& jsonResponse, std::string const& defaultMsg = "", uint16_type statusCode = invalid_uint16_type_value );

      private:
        std::shared_ptr<WorldComm> M_worldComm;
        std::string M_url, M_apiKey, M_token, M_path;
        std::set<std::string> M_fileIds, M_folderIds, M_itemIds;
    };

    struct FolderInfo : public std::tuple<std::string, std::string, size_type>
    {
        using super = std::tuple<std::string, std::string, size_type>;
        FolderInfo( std::string const& name = "", std::string const& id = "", size_type size = invalid_v<size_type> )
            : super( name, id, size )
        {
        }
        FolderInfo( FolderInfo const& ) = default;
        FolderInfo( FolderInfo&& ) = default;
        FolderInfo& operator=( FolderInfo const& ) = default;
        FolderInfo& operator=( FolderInfo&& ) = default;

        std::string const& name() const { return std::get<0>( *this ); }
        std::string const& id() const { return std::get<1>( *this ); }
        size_type size() const { return std::get<2>( *this ); }

        void addFolder( std::shared_ptr<FolderInfo> const& folder ) { M_folders.push_back( folder ); }
        void addItem( std::shared_ptr<ItemInfo> const& item ) { M_items.push_back( item ); }
        void addFile( std::shared_ptr<FileInfo> const& file ) { M_files.push_back( file ); }

        std::ostringstream print( size_t nTab = 0 ) const;

      private:
        std::vector<std::shared_ptr<FolderInfo>> M_folders;
        std::vector<std::shared_ptr<ItemInfo>> M_items;
        std::vector<std::shared_ptr<FileInfo>> M_files;
    };

    struct ItemInfo : public std::tuple<std::string, std::string, size_type>
    {
        using super = std::tuple<std::string, std::string, size_type>;
        ItemInfo( std::string const& name = "", std::string const& id = "", size_type size = invalid_v<size_type> )
            : super( name, id, size )
        {
        }
        ItemInfo( ItemInfo const& ) = default;
        ItemInfo( ItemInfo&& ) = default;
        ItemInfo& operator=( ItemInfo const& ) = default;
        ItemInfo& operator=( ItemInfo&& ) = default;

        std::string const& name() const { return std::get<0>( *this ); }
        std::string const& id() const { return std::get<1>( *this ); }
        size_type size() const { return std::get<2>( *this ); }

        void add( std::shared_ptr<FileInfo> const& file ) { M_files.push_back( file ); }

        std::ostringstream print( size_t nTab = 0 ) const;

      private:
        std::vector<std::shared_ptr<FileInfo>> M_files;
    };

    struct FileInfo : public std::tuple<std::string, std::string, size_type>
    {
        using super = std::tuple<std::string, std::string, size_type>;
        FileInfo( std::string const& name = "", std::string const& id = "", size_type size = invalid_v<size_type> )
            : super( name, id, size )
        {
        }
        FileInfo( FileInfo const& ) = default;
        FileInfo( FileInfo&& ) = default;
        FileInfo& operator=( FileInfo const& ) = default;
        FileInfo& operator=( FileInfo&& ) = default;

        std::string const& name() const { return std::get<0>( *this ); }
        std::string const& id() const { return std::get<1>( *this ); }
        size_type size() const { return std::get<2>( *this ); }
        std::string const& mimeType() const { return M_mimeType; }
        std::string const& checksum() const { return M_checksum; }
        std::string const& checksumType() const { return M_checksumType; }

        std::ostringstream print( size_t nTab = 0, bool isFileInItem = false ) const;

        void setMimeType( std::string const& s ) { M_mimeType = s; }
        void setChecksum( std::string const& type, std::string const& value )
        {
            M_checksumType = type;
            M_checksum = value;
        }

      private:
        std::string M_mimeType, M_checksumType, M_checksum;
    };

  private:
    std::shared_ptr<WorldComm> M_worldComm;
    boost::optional<URL> M_url;
    boost::optional<Github> M_github;
    boost::optional<Girder> M_girder;
};

} // namespace Feel

#endif