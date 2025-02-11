//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
//! @date 2024-12-01
//! @copyright 2024 Feel++ Consortium
//!

#include <boost/algorithm/string.hpp>
#include <cpr/cpr.h>
#include <feel/feelcore/remotedata.hpp>
#include <fstream>
#include <regex>

namespace Feel
{
RemoteData::Girder::Girder( std::string const& desc, WorldComm& worldComm )
    : M_worldComm( worldComm.shared_from_this() )
{
    std::regex ex("([ ]*)girder([ ]*):([ ]*)([{])([^]*)([}])");
    std::cmatch what;
    if( !regex_match(desc.c_str(), what, ex) )
        return;

    CHECK( what.size() == 7 ) << "invalid size";
    std::vector<std::string> keysvalues;
    std::string exprtosplit = std::string(what[5].first, what[5].second);
    //std::cout << fmt::format( "exprtosplit: {}", exprtosplit ) << "\n";
    auto resConvertion = convertDescToJson( exprtosplit );
    if ( !resConvertion.first )
    {
        if ( M_worldComm->isMasterRank() )
            std::cout << "Girder desc has a syntax error: " << desc << "\n";
        return;
    }
    nl::json jsonObj = resConvertion.second;

    if ( jsonObj.contains( "url" ) )
        M_url = jsonObj["url"].get<std::string>();
    if ( jsonObj.contains( "api_key" ) )
        M_apiKey = jsonObj["api_key"].get<std::string>();
    if ( jsonObj.contains( "token" ) )
        M_token = jsonObj["token"].get<std::string>();

    if ( jsonObj.contains( "file" ) )
    {
        if ( jsonObj["file"].is_array() )
        {
            for ( const auto& item : jsonObj["file"] )
                M_fileIds.insert( item.get<std::string>() );
        }
        else
        {
            M_fileIds.insert( jsonObj["file"].get<std::string>() );
        }
    }

    if ( jsonObj.contains( "folder" ) )
    {
        if ( jsonObj["folder"].is_array() )
        {
            for ( const auto& item : jsonObj["folder"] )
                M_folderIds.insert( item.get<std::string>() );
        }
        else
        {
            M_folderIds.insert( jsonObj["folder"].get<std::string>() );
        }
    }

    if ( jsonObj.contains( "item" ) )
    {
        if ( jsonObj["item"].is_array() )
        {
            for ( const auto& item : jsonObj["item"] )
                M_itemIds.insert( item.get<std::string>() );
        }
        else
        {
            M_itemIds.insert( jsonObj["item"].get<std::string>() );
        }
    }

    if ( jsonObj.contains( "path") )
    {
        //std::cout << fmt::format( "path: {}\n", jsonObj["path"].get<std::string>() );
        M_path = jsonObj["path"].get<std::string>();
    }

    if ( M_url.empty() )
        M_url = "https://girder.math.unistra.fr";

    if ( M_apiKey.empty() )
    {
        char* env = getenv( "FEELPP_GIRDER_API_KEY" );
        if ( env != NULL && env[0] != '\0' )
        {
            M_apiKey = env;
        }
    }

    if ( M_token.empty() )
    {
        char* env = getenv( "FEELPP_GIRDER_TOKEN" );
        if ( env != NULL && env[0] != '\0' )
        {
            M_token = env;
        }
    }

#if 0
    std::cout << "url: " << M_url << "\n";
    for (std::string const& fileId : M_fileIds)
        std::cout << "file id: " << fileId << "\n";
#endif
}

void RemoteData::Girder::setFolderIds( std::string const& folderId )
{
    M_folderIds.clear();
    M_folderIds.insert( folderId );
}

bool RemoteData::Girder::isInit() const
{
    return !M_url.empty();
}
bool RemoteData::Girder::canDownload() const
{
    return this->isInit() && ( !M_fileIds.empty() || !M_folderIds.empty() || !M_itemIds.empty() || !M_path.empty() );
}
bool RemoteData::Girder::canUpload() const
{
    return this->isInit() && ( !M_token.empty() || !M_apiKey.empty() );
}

nl::json RemoteData::Girder::getResourceInfoById(const std::string& resourceId, const std::string& token) const
{
    // List of possible resource types to check
    const std::vector<std::string> resourceTypes = { "item", "file", "folder" };

    for (const auto& resourceType : resourceTypes)
    {
        try
        {
            // Construct the URL for the resource endpoint with query parameters
            std::string url = fmt::format("{}/api/v1/resource/{}", M_url, resourceId);

            cpr::Response res = cpr::Get(
                cpr::Url{url},
                cpr::Header{ {"Girder-Token", token} },
                cpr::Parameters{
                    {"id", resourceId},
                    {"type", resourceType}
                },
                cpr::VerifySsl{false} // Disable SSL verification for debugging (enable in production)
            );

            // If the request succeeds, return the parsed JSON
            if (res.status_code == 200)
            {
                auto j = nl::json::parse(res.text);
                j["_modelType"] = resourceType;
                std::cout << fmt::format("Resource found with type '{}': {}", resourceType, j.dump()) << "\n";
                return j;
            }
        }
        catch (const std::exception& e)
        {
            // Log the exception for debugging but continue with the next resource type
            std::cerr << fmt::format("Error fetching resource with type '{}': {}", resourceType, e.what()) << "\n";
        }
    }

    // If all attempts fail, throw an exception
    throw std::runtime_error(fmt::format("Failed to find resource with ID: {}", resourceId));
}

std::vector<std::string>
RemoteData::Girder::download( const std::string& dir, const std::string& path ) const
{
    std::string token = (M_token.empty() && !M_apiKey.empty()) ? createToken() : std::string{};
    nl::json resourceInfo;
    try
    {
        resourceInfo = resourceLookup( path, token );
    }
    catch ( const std::exception& e )
    {
        std::cout << "Error in resourceLookup: " << e.what() << "\n";
        return {};
    }

    if ( !resourceInfo.contains( "_modelType" ) || !resourceInfo.contains( "_id" ) )
    {
        std::cout << "Invalid resource info: missing _modelType or _id\n";
        return {};
    }

    std::string resourceType = resourceInfo["_modelType"].get<std::string>();
    std::string resourceId = resourceInfo["_id"].get<std::string>();
    std::string name = resourceInfo["name"];
    //std::cout << fmt::format( "Downloading resource {} type: {} id: {}\n", name, resourceType, resourceId );
    if ( resourceType == "file" )
    {
        return { downloadFile( resourceId, dir, token ) };
    }
    else if ( resourceType == "folder" )
    {
        return { downloadFolder( resourceId, dir, token ) };
    }
    else if ( resourceType == "item" )
    {
        std::string url =  fmt::format("{}/api/v1/{}/{}/download",M_url,resourceType,resourceId);
        std::cout << fmt::format( "Downloading item {} from url: {}\n", name, url );
        cpr::Header headers;
        if ( !token.empty() )
            headers.insert( { "Girder-Token", token } );
        cpr::Response fileRes = cpr::Get(
            cpr::Url{url},
            headers,
            cpr::VerifySsl{false} // Add this if SSL verification causes issues
        );

        if (fileRes.status_code != 200) 
        {
            throw std::runtime_error(fmt::format("Failed to download file. HTTP status code: {}\n error message: {} ", fileRes.status_code, fileRes.text) );
        }

        // Save the file content to a local file
        std::string outputFileName = fmt::format("{}.zip",name); // Set desired output file name
        fs::path outputFilePath = fs::path(dir) / outputFileName;
        std::ofstream outputFile(outputFilePath, std::ios::binary);
        outputFile << fileRes.text;
        outputFile.close();

        //std::cout << "File downloaded successfully: " << outputFilePath << std::endl;
        return {outputFilePath.string()}; // Placeholder for item download implementation
    }
    else
    {
        throw std::runtime_error( "Unsupported resource type for download" );
    }
    return {};
}
std::vector<std::string>
RemoteData::Girder::download( std::string const& dir ) const
{
    std::vector<std::string> downloadedFileOrFolder;
    if ( M_worldComm->isMasterRank() )
    {
        if ( !fs::exists( dir ) )
            fs::create_directories( dir );
        // use token if given else create token if api key given
        std::string token = M_token;
        if ( M_token.empty() && !M_apiKey.empty() )
            token = this->createToken();
        // download girder files
        for ( std::string const& fileId : M_fileIds )
        {
            std::string file = this->downloadFile( fileId, dir, token );
            if ( !file.empty() )
                downloadedFileOrFolder.push_back( file );
        }
        // download girder folders
        for ( std::string const& folderId : M_folderIds )
        {
            std::string file = this->downloadFolder( folderId, dir, token );
            if ( !file.empty() )
                downloadedFileOrFolder.push_back( file );
        }
        // download girder items
        for ( std::string const& itemId : M_itemIds )
        {
            std::string file = this->downloadItem( itemId, dir, token );
            if ( !file.empty() )
                downloadedFileOrFolder.push_back( file );
        }
        if ( !M_path.empty() )
        {
            auto vs = download( dir, M_path );
            downloadedFileOrFolder.insert( downloadedFileOrFolder.end(), vs.begin(), vs.end() );
        }
        // delete token if created
        if ( M_token.empty() && !M_apiKey.empty() && !token.empty() )
            this->removeToken( token );
    }
    mpi::broadcast( M_worldComm->globalComm(), downloadedFileOrFolder, M_worldComm->masterRank() );
    return downloadedFileOrFolder;
}

std::string
RemoteData::Girder::errorMessage( nl::json const& jsonResponse, std::string const& defaultMsg, uint16_type statusCode )
{
    std::string errMsg = defaultMsg;
    if ( jsonResponse.contains( "message" ) )
        errMsg = jsonResponse["message"].get<std::string>();
    std::string errType;
    if ( jsonResponse.contains( "type" ) )
        errType = jsonResponse["type"].get<std::string>();
    std::string errMsgInfo;
    if ( statusCode != invalid_uint16_type_value )
        errMsgInfo += ( boost::format( "code %1%" ) % statusCode ).str();
    if ( !errMsgInfo.empty() )
        errMsgInfo += ",";
    if ( !errType.empty() )
        errMsgInfo += ( boost::format( "type %1%" ) % errType ).str();
    std::string res = "Girder error";
    if ( !errMsgInfo.empty() )
        res += " (" + errMsgInfo + ")";
    if ( !errMsg.empty() )
        res += " : " + errMsg;
    return res;
}

std::string
RemoteData::Girder::downloadFile( std::string const& fileId, std::string const& dir, std::string const& token ) const
{
    std::string downloadedFile;
    // Get metadata info
    std::string urlFileInfo = M_url + "/api/v1/file/" + fileId;
    std::vector<std::string> headersFileInfo;
    headersFileInfo.push_back( "Accept: application/json" );
    if ( !token.empty() )
        headersFileInfo.push_back( "Girder-Token: " + token );
    std::ostringstream omemfile;
    StatusRequestHTTP status = requestHTTPGET( urlFileInfo, headersFileInfo, omemfile );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPGET : " << status.msg() << "\n";
        return {};
    }
    // Parse JSON
    nl::json jsonResponse = nl::json::parse( omemfile.str() );
    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( jsonResponse, "Getting metadata (before download) fails", status.code() ) << "\n";
        return {};
    }

    // Extract info from JSON
    if ( !jsonResponse.contains( "name" ) )
    {
        std::cerr << "Invalid ID: Not a file or does not exist\n";
        return {};
    }
    std::string filename = jsonResponse["name"].get<std::string>();
    std::string mimeType = jsonResponse.value( "mimeType", "" );
    std::string sha512 = jsonResponse.value( "sha512", "" );

    std::string filepath = ( fs::path( dir ) / filename ).string();
    std::string metadatapath = ( fs::path( dir ) / ( filename + ".metadata.json" ) ).string();

    // Check if download is necessary
    bool doDownload = true;
    if ( !sha512.empty() && fs::exists( filepath ) && fs::is_regular_file( filepath ) && fs::exists( metadatapath ) )
    {
        nl::json existingMetadata;
        std::ifstream metadataFile( metadatapath );
        metadataFile >> existingMetadata;
        if ( existingMetadata.contains( "sha512" ) )
        {
            if ( sha512 == existingMetadata["sha512"].get<std::string>() )
                doDownload = false;
        }
    }
    // Download the file
    if ( doDownload )
    {
        std::string urlFileDownload = M_url + "/api/v1/file/" + fileId + "/download";
        std::vector<std::string> headersFileDownload;
        if ( !mimeType.empty() )
            headersFileDownload.push_back( "Accept: " + mimeType );
        if ( !token.empty() )
            headersFileDownload.push_back( "Girder-Token: " + token );
        std::ofstream ofile( filepath, std::ios::out | std::ios::binary );
        status = requestHTTPGET( urlFileDownload, headersFileDownload, ofile );
        ofile.close();
        if ( !status.success() )
        {
            std::cout << "Girder error in requestHTTPGET : " << status.msg() << "\n";
            return {};
        }
        if ( status.code() != 200 )
        {
            std::cout << Girder::errorMessage( nl::json{}, "Downloading file fails", status.code() ) << "\n";
            return {};
        }
        // Save metadata
        std::ofstream ofileMetadata( metadatapath, std::ios::out );
        ofileMetadata << omemfile.str();
        ofileMetadata.close();
    }

    downloadedFile = filepath;
    return downloadedFile;
}
std::string
RemoteData::Girder::downloadFolder( std::string const& folderId, std::string const& dir, std::string const& token ) const
{
    std::string downloadedFolder;
    // Get metadata info
    std::string urlFolderInfo = M_url + "/api/v1/folder/" + folderId;
    std::vector<std::string> headersFolderInfo;
    headersFolderInfo.push_back( "Accept: application/json" );
    if ( !token.empty() )
        headersFolderInfo.push_back( "Girder-Token: " + token );
    std::ostringstream omemfile;
    StatusRequestHTTP status = requestHTTPGET( urlFolderInfo, headersFolderInfo, omemfile );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPGET: " << status.msg() << "\n";
        return {};
    }

    // Parse the JSON response
    nl::json jsonResponse;
    try
    {
        jsonResponse = nl::json::parse( omemfile.str() );
    }
    catch ( nl::json::parse_error& e )
    {
        std::cout << "Error parsing JSON response: " << e.what() << "\n";
        return {};
    }

    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( jsonResponse, "Getting metadata (before download) fails", status.code() ) << "\n";
        return {};
    }

    // Extract folder name from jsonResponse
    if ( !jsonResponse.contains( "name" ) )
    {
        std::cout << "Invalid ID: Not a folder or does not exist\n";
        return {};
    }
    std::string foldername = jsonResponse["name"].get<std::string>();

    // Download the folder
    std::string urlFolderDownload = M_url + "/api/v1/folder/" + folderId + "/download";
    std::vector<std::string> headersFolderDownload;
    if ( !token.empty() )
        headersFolderDownload.push_back( "Girder-Token: " + token );
    std::string filepath = ( fs::path( dir ) / ( foldername + ".zip" ) ).string();
    std::ofstream ofile( filepath, std::ios::out | std::ios::binary );
    status = requestHTTPGET( urlFolderDownload, headersFolderDownload, ofile );
    ofile.close();
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPGET: " << status.msg() << "\n";
        return {};
    }
    if ( status.code() != 200 )
    {
        // Since we don't have a JSON response for the download error, pass an empty json object
        std::cout << Girder::errorMessage( nl::json{}, "Downloading folder fails", status.code() ) << "\n";
        return {};
    }

    // Save metadata
    std::string metadatapath = ( fs::path( dir ) / ( foldername + ".metadata.json" ) ).string();
    std::ofstream ofileMetadata( metadatapath, std::ios::out );
    ofileMetadata << jsonResponse.dump( 4 ); // Write the JSON with indentation
    ofileMetadata.close();

    downloadedFolder = filepath;
    return downloadedFolder;
}

std::string
RemoteData::Girder::downloadItem(const std::string& resourceId, const std::string& dir, const std::string& token) const
{
    // Construct the URL to get item metadata
    std::string urlItemInfo = M_url + "/api/v1/item/" + resourceId;
    std::vector<std::string> headersItemInfo;
    headersItemInfo.push_back("Accept: application/json");
    if (!token.empty())
        headersItemInfo.push_back("Girder-Token: " + token);
    
    // Fetch item metadata
    std::ostringstream omemfile;
    StatusRequestHTTP status = requestHTTPGET(urlItemInfo, headersItemInfo, omemfile);
    if (!status.success())
    {
        std::cout << "Girder error in requestHTTPGET: " << status.msg() << "\n";
        return {};
    }

    // Parse the JSON response
    nl::json jsonResponse;
    try
    {
        jsonResponse = nl::json::parse(omemfile.str());
    }
    catch (nl::json::parse_error& e)
    {
        std::cout << "Error parsing JSON response: " << e.what() << "\n";
        return {};
    }

    if (status.code() != 200)
    {
        std::cout << Girder::errorMessage(jsonResponse, "Getting metadata (before download) fails", status.code()) << "\n";
        return {};
    }

    // Extract item name from metadata
    if (!jsonResponse.contains("name"))
    {
        std::cout << "Invalid ID: Not an item or does not exist\n";
        return {};
    }
    std::string itemName = jsonResponse["name"].get<std::string>();

    // Construct the URL to download the item
    std::string urlItemDownload = M_url + "/api/v1/item/" + resourceId + "/download";
    std::vector<std::string> headersItemDownload;
    if (!token.empty())
        headersItemDownload.push_back("Girder-Token: " + token);
    
    // Download the item
    std::string filepath = (fs::path(dir) / (itemName + ".zip")).string();
    std::ofstream ofile(filepath, std::ios::out | std::ios::binary);
    status = requestHTTPGET(urlItemDownload, headersItemDownload, ofile);
    ofile.close();
    if (!status.success())
    {
        std::cout << "Girder error in requestHTTPGET: " << status.msg() << "\n";
        return {};
    }
    if (status.code() != 200)
    {
        std::cout << Girder::errorMessage(nl::json{}, "Downloading item fails", status.code()) << "\n";
        return {};
    }

    // Save item metadata
    std::string metadatapath = (fs::path(dir) / (itemName + ".metadata.json")).string();
    std::ofstream ofileMetadata(metadatapath, std::ios::out);
    ofileMetadata << jsonResponse.dump(4); // Write the JSON with indentation
    ofileMetadata.close();

    // Return the downloaded file path
    return filepath;
}

std::vector<std::string>
RemoteData::Girder::upload(const std::string& dataPath, const std::string& destination, bool sync) const
{
    auto res = this->upload(std::vector<std::pair<std::string, std::string>>(1, std::make_pair(dataPath, destination)), sync);
    CHECK(res.size() == 1) << "Wrong size " << res.size() << " : must be 1";
    return res[0];
}

std::vector<std::vector<std::string>>
RemoteData::Girder::upload(const std::vector<std::pair<std::string, std::string>>& dataToUpload, bool sync) const
{
    if (dataToUpload.empty())
        return {};

    CHECK(canUpload()) << "Authentication unavailable";

    std::vector<std::vector<std::string>> res;
    res.reserve(dataToUpload.size());

    if (M_worldComm->isMasterRank())
    {
        // Use token if given else create token if API key given
        std::string token = M_token;
        if (M_token.empty())
            token = this->createToken();

        for ( auto [dataPath,destination] : dataToUpload )
        {
            if (!fs::exists(dataPath))
            {
                std::cout << "Warning in Girder upload, data path does not exist: " << dataPath << "\n";
                continue;
            }
            std::string resourceType;
            std::string resourceId;
            std::cout << fmt::format("Uploading data from {} to {}\n", dataPath, destination) << "\n";
            // Check if destination starts with '/', indicating a path
            if (!destination.empty() && destination[0] == '/')
            {
                // It's a path, use resourceLookup
                nl::json resourceInfo;
                try
                {
                    resourceInfo = resourceLookup(destination, token);
                    std::cout << fmt::format("resourceInfo: {}\n", resourceInfo.dump(4)) << "\n";
                }
                catch (const std::exception& e)
                {
                    std::cout << "Error in resourceLookup: " << e.what() << "\n";
                    continue;
                }

                if (!resourceInfo.contains("_modelType") || !resourceInfo.contains("_id"))
                {
                    std::cout << "Invalid resource info: missing _modelType or _id\n";
                    continue;
                }

                resourceType = resourceInfo["_modelType"].get<std::string>();
                resourceId = resourceInfo["_id"].get<std::string>();
            }
            else
            {
                // Destination is an ID; determine the resource type
                try
                {
                    nl::json resourceInfo = getResourceInfoById(destination, token);
                    if (!resourceInfo.contains("_modelType") || !resourceInfo.contains("_id"))
                    {
                        std::cout << "Invalid resource info: missing _modelType or _id\n";
                        continue;
                    }

                    resourceType = resourceInfo["_modelType"].get<std::string>();
                    resourceId = resourceInfo["_id"].get<std::string>();
                }
                catch (const std::exception& e)
                {
                    std::cout << "Error fetching resource info by ID: " << e.what() << "\n";
                    continue;
                }
            }

            fs::path dataFsPath(dataPath);
            if (resourceType == "folder")
            {
                // Upload to folder
                std::vector<std::string> uploadedResources;
                uploadDirectoryToFolder(dataPath, resourceId, token, uploadedResources);
                res.push_back( uploadedResources );
            }
            else if (resourceType == "item")
            {
                // Upload to item
                std::vector<std::string> uploadedResources;
                uploadFilesToItem(dataPath, resourceId, token, uploadedResources);
                res.push_back( uploadedResources );
            }
            else
            {
                std::cout << "Unsupported resource type for upload: " << resourceType << "\n";
                continue;
            }
        }

        // Delete token if created
        if (M_token.empty() && !token.empty())
            this->removeToken(token);
    }

    if (sync)
        mpi::broadcast(M_worldComm->globalComm(), res, M_worldComm->masterRank());

    return res;
}

std::string
RemoteData::Girder::createItemImpl(const std::string& itemName, const std::string& parentFolderId, const std::string& token) const
{
    // Construct the URL for creating an item
    std::string urlCreateItem = M_url + "/api/v1/item";
    std::cout << fmt::format("Creating item {} in folder {} using token: {}\n", itemName, parentFolderId, token) << "\n";
    // Set up the headers
    cpr::Header headers = {
       // {"Accept", "application/json"},
       // {"Content-Type", "application/json"},
        {"Girder-Token", token}
    };
    
    // Send the POST request to create the item
    cpr::Response res = cpr::Post(
        cpr::Url{urlCreateItem},
        headers,
        cpr::Payload{
            {"name", itemName},
            {"folderId", parentFolderId},
            {"reuseExisting", std::to_string(true)}
        },
        cpr::VerifySsl{false} // Set to true in production
    );
    std::cout << fmt::format("Create item response: {}\n", res.text) << "\n";
    // Check if the request was successful
    if (res.status_code != 200)
    {
        // Parse the JSON response to extract error messages
        nl::json jsonResponse;
        try
        {
            jsonResponse = nl::json::parse(res.text);
        }
        catch (nl::json::parse_error& e)
        {
            std::cout << "Error parsing JSON response: " << e.what() << "\n";
            return {};
        }
        // Output the error message
        std::cout << Girder::errorMessage(jsonResponse, "Create item fails", res.status_code) << "\n";
        return {};
    }
    
    // Parse the successful response to get the item ID
    nl::json jsonResponse;
    try
    {
        jsonResponse = nl::json::parse(res.text);
    }
    catch (nl::json::parse_error& e)
    {
        std::cout << "Error parsing JSON response: " << e.what() << "\n";
        return {};
    }
    
    // Extract the item ID from the response
    std::string itemId = jsonResponse.value("_id", "");
    return itemId;
}


std::vector<std::pair<std::string, std::string>>
RemoteData::Girder::createItem( std::string const& itemPath, std::string const& parentId, bool sync ) const
{
    
    std::vector<std::tuple<std::string, std::string>> itemInfo;
    if ( M_worldComm->isMasterRank() )
    {
        std::string token = ( M_token.empty() && !M_apiKey.empty() )
                                 ? createToken()
                                 : std::string{};
        std::cout << fmt::format( "Creating item from {} in folder {} using toekn: {}\n", itemPath, parentId, token ) << "\n";
        std::string itemId = createItemImpl( fs::path( itemPath ).filename().string(), parentId, token );
        if ( !itemId.empty() )
            itemInfo.push_back( std::make_tuple( itemPath, itemId ) );
    }
    if ( sync )
        mpi::broadcast( M_worldComm->globalComm(), itemInfo, M_worldComm->masterRank() );
    std::vector<std::pair<std::string, std::string>> res;
    for ( auto [itemPath, itemId] : itemInfo )
        res.push_back( std::make_pair( itemPath, itemId ) );
    return res;
}

std::string
RemoteData::Girder::createFolderImpl(const std::string& folderName, const std::string& parentId, const std::string& token, const std::string& parentType ) const
{
    // Construct the URL for creating a folder
    std::string urlCreateFolder = M_url + "/api/v1/folder";
    
    // Set up the headers
    cpr::Header headers = {
        {"Accept", "application/json"},
        {"Content-Type", "application/json"},
        {"Girder-Token", token}
    };
    
    // Create the JSON body for the POST request
    nl::json requestBody = {
        {"name", folderName},
        {"parentId", parentId},
        {"parentType", parentType},
        {"reuseExisting", true}
    };
    
    // Send the POST request to create the folder
    cpr::Response res = cpr::Post(
        cpr::Url{urlCreateFolder},
        headers,
        cpr::Body{requestBody.dump()},
        cpr::VerifySsl{false} // Set to true in production
    );
    
    // Check if the request was successful
    if (res.status_code != 200)
    {
        // Parse the JSON response to extract error messages
        nl::json jsonResponse;
        try
        {
            jsonResponse = nl::json::parse(res.text);
        }
        catch (nl::json::parse_error& e)
        {
            std::cout << "Error parsing JSON response: " << e.what() << "\n";
            return {};
        }
        // Output the error message
        std::cout << Girder::errorMessage(jsonResponse, "Create folder fails", res.status_code) << "\n";
        return {};
    }
    
    // Parse the successful response to get the folder ID
    nl::json jsonResponse;
    try
    {
        jsonResponse = nl::json::parse(res.text);
    }
    catch (nl::json::parse_error& e)
    {
        std::cout << "Error parsing JSON response: " << e.what() << "\n";
        return {};
    }
    
    // Extract the folder ID from the response
    std::string folderId = jsonResponse.value("_id", "");
    return folderId;
}
void
RemoteData::Girder::uploadDirectoryToFolder(const std::string& localDir, const std::string& parentFolderId, const std::string& token, std::vector<std::string>& uploadedResources) const
{
    fs::path dataFsPath(localDir);
    if (fs::is_regular_file(dataFsPath))
    {
        // Upload the file to a new item within the folder
        std::string itemName = dataFsPath.filename().string();

        std::string itemId = createItemImpl(itemName, parentFolderId, token);
        if (!itemId.empty())
        {
            std::string fileId = uploadFileImpl(dataFsPath.string(), itemId, token, "item");
            if (!fileId.empty())
            {
                uploadedResources.push_back(fileId);
            }
        }
    }
    else if (fs::is_directory(dataFsPath))
    {
        // Create a folder in Girder
        std::string folderName = dataFsPath.filename().string();
        std::string subfolderId = createFolderImpl(folderName, parentFolderId, token);
        if (!subfolderId.empty())
        {
            // Recursively upload the contents of the subfolder
            for (const fs::directory_entry& entry : fs::directory_iterator(dataFsPath))
            {
                uploadDirectoryToFolder(entry.path().string(), subfolderId, token, uploadedResources);
            }
        }
    }
    else
    {
        std::cout << "Unsupported file system object: " << localDir << "\n";
    }
}
void
RemoteData::Girder::uploadFilesToItem(const std::string& localPath, const std::string& itemId, const std::string& token, std::vector<std::string>& uploadedResources) const
{
    fs::path dataFsPath(localPath);
    if (fs::is_regular_file(dataFsPath))
    {
        // Upload the file to the item
        std::string fileId = uploadFileImpl(dataFsPath.string(), itemId, token, "item");
        if (!fileId.empty())
        {
            uploadedResources.push_back(fileId);
        }
    }
    else if (fs::is_directory(dataFsPath))
    {
        // Optionally handle subdirectories
        for (const fs::directory_entry& entry : fs::directory_iterator(dataFsPath))
        {
            uploadFilesToItem(entry.path().string(), itemId, token, uploadedResources);
        }
    }
    else
    {
        std::cout << "Unsupported file system object: " << localPath << "\n";
    }
}
void RemoteData::Girder::replaceFile( std::string const& filePath, std::string const& fileId ) const
{
    this->replaceFile( std::vector<std::pair<std::string, std::string>>( 1, std::make_pair( filePath, fileId ) ) );
}

void RemoteData::Girder::replaceFile( std::vector<std::pair<std::string, std::string>> const& filesToReplace ) const
{
    if ( filesToReplace.empty() )
        return;

    CHECK( !M_token.empty() || !M_apiKey.empty() ) << "authentication unavailable";

    if ( M_worldComm->isMasterRank() )
    {
        // use token if given else create token if api key given
        std::string token = M_token;
        if ( M_token.empty() )
            token = this->createToken();
        // call replace file request for each file
        for ( int k = 0; k < filesToReplace.size(); ++k )
        {
            std::string filePath = filesToReplace[k].first;
            std::string fileId = filesToReplace[k].second;
            this->replaceFileImpl( filePath, fileId, token );
        }
        // delete token if created
        if ( M_token.empty() && !token.empty() )
            this->removeToken( token );
    }
}

std::vector<std::pair<std::string, std::string>>
RemoteData::Girder::createFolder( std::string const& folderPath, std::string const& parentId, bool sync ) const
{
    CHECK( !M_token.empty() || !M_apiKey.empty() ) << "authentication unavailable";
    std::string currentParentId = parentId;
    if ( parentId.empty() && !M_folderIds.empty() )
        currentParentId = *M_folderIds.begin();
    CHECK( !currentParentId.empty() ) << "a parentId is required";

    std::vector<boost::tuple<std::string, std::string>> foldersInfo;
    if ( M_worldComm->isMasterRank() )
    {
        std::string token = M_token;
        if ( M_token.empty() )
            token = this->createToken();

        fs::path folderFsPath( folderPath );
        for ( auto const& subdir : folderFsPath )
        {
            if ( !Feel::filename_is_dot( subdir ) )
            {
                currentParentId = this->createFolderImpl( subdir.string(), currentParentId, token );
                foldersInfo.push_back( boost::make_tuple( subdir.string(), currentParentId ) );
            }
        }

        if ( M_token.empty() && !token.empty() )
            this->removeToken( token );
    }
    if ( sync )
        mpi::broadcast( M_worldComm->globalComm(), foldersInfo, M_worldComm->masterRank() );
    std::vector<std::pair<std::string, std::string>> res;
    for ( auto const& folderInfo : foldersInfo )
        res.push_back( std::make_pair( boost::get<0>( folderInfo ), boost::get<1>( folderInfo ) ) );
    return res;
}

std::vector<std::string>
RemoteData::Girder::uploadRecursively( std::string const& dataPath, std::string const& parentId, std::string const& token ) const
{
    std::vector<std::string> res;
    fs::path dataFsPath( dataPath );
    if ( fs::is_directory( dataFsPath ) )
    {
        std::string folderName = dataFsPath.filename().string();
        std::string newParentId = this->createFolderImpl( folderName, parentId, token );

        fs::directory_iterator end_itr;
        for ( fs::directory_iterator itr( dataFsPath ); itr != end_itr; ++itr )
        {
            auto filesUploaded = this->uploadRecursively( itr->path().string(), newParentId, token );
            for ( std::string const& fileUploaded : filesUploaded )
                if ( !fileUploaded.empty() )
                    res.push_back( fileUploaded );
        }
    }
    else if ( fs::is_regular_file( dataFsPath ) )
    {
        std::string fileUploaded = this->uploadFileImpl( dataPath, parentId, token, "folder" );
        if ( !fileUploaded.empty() )
            res.push_back( fileUploaded );
    }
    return res;
}
#if 1
std::string
RemoteData::Girder::uploadFileImpl(const std::string& filepath, const std::string& parentId, const std::string& token, const std::string& parentType) const
{
    // Ensure parentType is either "folder" or "item"
    if (parentType != "folder" && parentType != "item")
    {
        throw std::invalid_argument("Invalid parentType. Must be 'folder' or 'item'.");
    }

    // Extract filename and file size
    fs::path filePathObj(filepath);
    std::string filename = filePathObj.filename().string();

    // Open the file
    std::ifstream fileStream(filepath, std::ios::binary);
    if (!fileStream)
    {
        std::cerr << "Failed to open file: " << filepath << "\n";
        return {};
    }

    // Get file size
    fileStream.seekg(0, std::ios::end);
    std::streamsize fileSize = fileStream.tellg();
    fileStream.seekg(0, std::ios::beg);

    // Initialize the upload
    std::string uploadId;
    try
    {
        uploadId = initializeUpload(filename, fileSize, parentId, parentType, token);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error initializing upload: " << e.what() << "\n";
        return {};
    }

    std::cout << fmt::format("Uploading file {} ({} bytes) to parent ID {} (type: {})\n", filename, fileSize, parentId, parentType);

    // Upload the file in chunks
    const std::streamsize chunkSize = 64 * 1024 * 1024; // 64 MB
    std::streamsize offset = 0;
    std::vector<char> buffer(chunkSize);
    std::string fileId;
    while (offset < fileSize)
    {
        std::streamsize bytesToRead = std::min(chunkSize, fileSize - offset);
        fileStream.read(buffer.data(), bytesToRead);

        try
        {
            std::cout << fmt::format("Uploading chunk at offset {} ({} bytes)\n", offset, bytesToRead);
            fileId = uploadChunk(uploadId, offset, buffer.data(), bytesToRead, token);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Error uploading chunk at offset " << offset << ": " << e.what() << "\n";
            return {};
        }

        offset += bytesToRead;
    }
    
    std::cout << "Upload complete\n";
    // After the last chunk, the response includes the file ID
    //std::string fileId = uploadChunk(uploadId, offset, nullptr, 0, token, true);
    std::cout << fmt::format("File uploaded successfully. File ID: {}\n", fileId);
    return fileId;
}
std::string
RemoteData::Girder::initializeUpload(const std::string& filename, std::streamsize fileSize, const std::string& parentId, const std::string& parentType, const std::string& token) const
{
    std::string url = M_url + "/api/v1/file";

    // Request body
    nl::json requestBody = {
        {"name", filename},
        {"parentId", parentId},
        {"parentType", parentType},
        {"size", fileSize}
    };
    cpr::Parameters params = {
        {"name", filename},
        {"parentId", parentId},
        {"parentType", parentType},
        {"size", std::to_string(fileSize)}
    };
    // Headers
    cpr::Header headers = {
        {"Accept", "application/json"},
        //{"Content-Type", "application/json"},
        {"Girder-Token", token}
    };

    // POST request to initialize upload
    cpr::Response res = cpr::Post(
        cpr::Url{url},
        params,
        headers,
        cpr::Body{requestBody.dump()},
        cpr::VerifySsl{false} // Set to true in production
    );

    if (res.status_code != 200)
    {
        throw std::runtime_error("Failed to initialize upload: " + res.text);
    }

    nl::json jsonResponse = nl::json::parse(res.text);
    std::string uploadId = jsonResponse.value("_id", "");

    if (uploadId.empty())
    {
        throw std::runtime_error("Upload ID not found in response.");
    }

    return uploadId;
}
std::string
RemoteData::Girder::uploadChunk(const std::string& uploadId, std::streamsize offset, const char* data, std::streamsize size, const std::string& token, bool isFinalChunk) const
{
    std::string url = M_url + "/api/v1/file/chunk";

    // Parameters
    cpr::Parameters params = {
        {"uploadId", uploadId},
        {"offset", std::to_string(offset)}
    };

    // Headers
    cpr::Header headers = {
        {"Accept", "application/json"},
        {"Girder-Token", token}
    };

    cpr::Response res;

    if (size > 0)
    {
        // Data to upload
        cpr::Body body(std::string(data, size));

        // Content-Type header
        headers["Content-Type"] = "application/octet-stream";

        // POST request to upload chunk
        res = cpr::Post(
            cpr::Url{url},
            params,
            headers,
            body,
            cpr::VerifySsl{false} // Set to true in production
        );
    }
    else if (isFinalChunk)
    {
        // Finalize the upload without data
        res = cpr::Post(
            cpr::Url{url},
            params,
            headers,
            cpr::VerifySsl{false} // Set to true in production
        );
    }
    else
    {
        throw std::invalid_argument("Invalid chunk upload: size is zero and not final chunk.");
    }

    if (res.status_code != 200)
    {
        throw std::runtime_error("Failed to upload chunk: " + res.text);
    }

    nl::json jsonResponse = nl::json::parse(res.text);
    std::cout << fmt::format("response: {}", jsonResponse.dump(4)) << "\n";
    if ( jsonResponse.contains( "_id" ) )
    {
        std::string fileId = jsonResponse["_id"].get<std::string>();
        return fileId;
    }
    // For intermediate chunks, return empty string
    return "";
}
#else
std::string
RemoteData::Girder::uploadFileImpl( std::string const& filepath, std::string const& parentId, std::string const& token ) const
{
    std::string filename = fs::path( filepath ).filename().string();
    std::string fileExtension = fs::path( filepath ).extension().string();

    std::ifstream imemfile( filepath, std::ios::binary );
    // Get length of file
    imemfile.seekg( 0, imemfile.end );
    long fsize = imemfile.tellg();
    imemfile.seekg( 0, imemfile.beg );

    std::string urlFileUpload = M_url + "/api/v1/file?parentType=folder&parentId=" + parentId;
    urlFileUpload += "&name=" + boost::replace_all_copy( filename, " ", "+" );
    urlFileUpload += "&size=" + std::to_string( fsize );
    std::string mimeType;
    auto itFindMimeType = mimeTypes.find( fileExtension );
    if ( itFindMimeType != mimeTypes.end() )
    {
        std::vector<std::string> exprSplitted;
        boost::split( exprSplitted, itFindMimeType->second, boost::is_any_of( "/" ), boost::token_compress_on );
        mimeType = exprSplitted[0];
        for ( size_t k = 1; k < exprSplitted.size(); ++k )
            mimeType += "%2F" + exprSplitted[k];
    }
    if ( !mimeType.empty() )
        urlFileUpload += "&mimeType=" + mimeType;

    std::vector<std::string> headersFileInfo;
    headersFileInfo.push_back( "Accept: application/json" );
    headersFileInfo.push_back( "Content-Type: application/form-data" );
    CHECK( !token.empty() ) << "A token is required for upload";
    headersFileInfo.push_back( "Girder-Token: " + token );

    long limitFileSizePosted = 67108864;
    long fsizePosted = std::min( fsize, limitFileSizePosted );
    headersFileInfo.push_back( "Content-Length: " + std::to_string( fsizePosted ) );

    std::ostringstream omemfile;
    StatusRequestHTTP status = requestHTTPPOST( urlFileUpload, headersFileInfo, imemfile, fsizePosted, omemfile );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPPOST: " << status.msg() << "\n";
        return {};
    }

    // Parse the response as JSON using nl::json
    nl::json jsonResponse;
    try
    {
        jsonResponse = nl::json::parse( omemfile.str() );
    }
    catch ( nl::json::parse_error& e )
    {
        std::cout << "Error parsing JSON response: " << e.what() << "\n";
        return {};
    }

    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( jsonResponse, "Upload file fails", status.code() ) << "\n";
        return {};
    }

    std::string fileIdCreated;
    if ( jsonResponse.contains( "_id" ) )
        fileIdCreated = jsonResponse["_id"].get<std::string>();

    long fileSizeStayToUpload = fsize - fsizePosted;
    if ( fileSizeStayToUpload > 0 )
    {
        std::vector<std::pair<long, long>> chunkDivisions;
        while ( fileSizeStayToUpload > 0 )
        {
            long fsizeToPostInChunk = std::min( fileSizeStayToUpload, limitFileSizePosted );
            chunkDivisions.push_back( std::make_pair( fsizePosted, fsizeToPostInChunk ) );
            fileSizeStayToUpload -= fsizeToPostInChunk;
            fsizePosted += fsizeToPostInChunk;
        }

        for ( const auto& chunkDivision : chunkDivisions )
        {
            long fileOffset = chunkDivision.first;
            std::string fileOffsetStr = std::to_string( fileOffset );
            long fsizeToPostInChunk = chunkDivision.second;
            std::string urlFileChunk = M_url + "/api/v1/file/chunk?uploadId=" + fileIdCreated + "&offset=" + fileOffsetStr;

            imemfile.seekg( fileOffset, imemfile.beg );

            std::vector<std::string> headersFileChunk;
            headersFileChunk.push_back( "Accept: application/json" );
            headersFileChunk.push_back( "Content-Type: application/form-data" );
            headersFileChunk.push_back( "Girder-Token: " + token );
            headersFileChunk.push_back( "Content-Length: " + std::to_string( fsizeToPostInChunk ) );

            std::ostringstream omemfileChunk;
            status = requestHTTPPOST( urlFileChunk, headersFileChunk, imemfile, fsizeToPostInChunk, omemfileChunk );
            if ( !status.success() )
            {
                std::cout << "Girder error in requestHTTPPOST: " << status.msg() << "\n";
                return {};
            }

            // Parse the chunk response as JSON
            nl::json jsonResponseChunk;
            try
            {
                jsonResponseChunk = nl::json::parse( omemfileChunk.str() );
            }
            catch ( nl::json::parse_error& e )
            {
                std::cout << "Error parsing JSON response: " << e.what() << "\n";
                return {};
            }

            if ( status.code() != 200 )
            {
                std::cout << Girder::errorMessage( jsonResponseChunk, "Upload file fails", status.code() ) << "\n";
                return {};
            }
        }
    }

    return fileIdCreated;
}
#endif // 0
void RemoteData::Girder::replaceFileImpl( std::string const& filePath, std::string const& fileId, std::string const& token ) const
{
    CHECK( !token.empty() ) << "A token is required for upload";
    CHECK( fs::is_regular_file( filePath ) ) << "Must be a file";

    std::ifstream imemfile( filePath, std::ios::binary );
    // Get length of file
    imemfile.seekg( 0, imemfile.end );
    long fsize = imemfile.tellg();
    imemfile.seekg( 0, imemfile.beg );

    // Request HTTP PUT: contents
    std::string urlFileContents = M_url + "/api/v1/file/" + fileId + "/contents";
    urlFileContents += "?size=" + std::to_string( fsize );
    std::vector<std::string> headersFileContents{
        "Accept: application/json",
        "Content-Type: application/json",
        "Content-Length: 0",
        "Girder-Token: " + token };

    std::ostringstream omemfileFileContents;
    StatusRequestHTTP status = requestHTTPCUSTOM( "PUT", urlFileContents, headersFileContents, omemfileFileContents );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPCUSTOM (PUT): " << status.msg() << "\n";
        return;
    }

    // Parse JSON response
    nl::json jsonResponse;
    try
    {
        jsonResponse = nl::json::parse( omemfileFileContents.str() );
    }
    catch ( nl::json::parse_error& e )
    {
        std::cout << "Error parsing JSON response: " << e.what() << "\n";
        return;
    }

    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( jsonResponse, "Replace file (contents) fails", status.code() ) << "\n";
        return;
    }

    std::string uploadId;
    if ( jsonResponse.contains( "_id" ) )
        uploadId = jsonResponse["_id"].get<std::string>();
    CHECK( !uploadId.empty() ) << "No _id in metadata returned";

    // Request HTTP POST: chunk
    std::string urlFileChunk = M_url + "/api/v1/file/chunk";
    urlFileChunk += "?uploadId=" + uploadId;
    urlFileChunk += "&offset=0";
    std::vector<std::string> headersFileChunk{
        "Accept: application/json",
        "Content-Type: application/form-data",
        "Girder-Token: " + token };

    std::ostringstream omemfile;
    status = requestHTTPPOST( urlFileChunk, headersFileChunk, imemfile, fsize, omemfile );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPPOST: " << status.msg() << "\n";
        return;
    }
    if ( status.code() != 200 )
    {
        // Parse JSON response
        nl::json jsonResponseChunk;
        try
        {
            jsonResponseChunk = nl::json::parse( omemfile.str() );
        }
        catch ( nl::json::parse_error& e )
        {
            std::cout << "Error parsing JSON response: " << e.what() << "\n";
            return;
        }

        std::cout << Girder::errorMessage( jsonResponseChunk, "Upload file (chunk) fails", status.code() ) << "\n";
        return;
    }
}

#if 0
std::string
RemoteData::Girder::createFolderImpl( std::string const& folderName, std::string const& parentId, std::string const& token ) const
{
    std::string urlCreateFolder = M_url + "/api/v1/folder?parentType=folder&parentId=" + parentId;
    urlCreateFolder += "&name=" + folderName;
    urlCreateFolder += "&reuseExisting=true";

    std::vector<std::string> headers{
        "Accept: application/json",
        "Content-Type: application/x-www-form-urlencoded",
        "Girder-Token: " + token };

    std::ostringstream omemfile;
    StatusRequestHTTP status = requestHTTPPOST( urlCreateFolder, headers, omemfile );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPPOST: " << status.msg() << "\n";
        return {};
    }

    // Parse JSON response
    nl::json jsonResponse;
    try
    {
        jsonResponse = nl::json::parse( omemfile.str() );
    }
    catch ( nl::json::parse_error& e )
    {
        std::cout << "Error parsing JSON response: " << e.what() << "\n";
        return {};
    }

    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( jsonResponse, "Create folder fails", status.code() ) << "\n";
        return {};
    }

    std::string folderIdCreated;
    if ( jsonResponse.contains( "_id" ) )
        folderIdCreated = jsonResponse["_id"].get<std::string>();
    return folderIdCreated;
}
#endif
std::string
RemoteData::Girder::createToken( int duration ) const
{
    std::string urlCreateToken = M_url + "/api/v1/api_key/token";
    urlCreateToken += "?key=" + M_apiKey;
    urlCreateToken += "&duration=" + std::to_string( duration );

    std::vector<std::string> headers{
        "Accept: application/json",
        "Content-Type: application/json" };

    std::ostringstream omemfile;
    StatusRequestHTTP status = requestHTTPPOST( urlCreateToken, headers, omemfile );
    if ( !status.success() )
    {
        std::cout << fmt::format("Girder error in createToken, requestHTTPPOST: {}",status.msg());
        throw std::runtime_error(fmt::format("Girder error in createToken, requestHTTPPOST: {}",status.msg()));
    }

    // Parse JSON response
    nl::json jsonResponse;
    try
    {
        jsonResponse = nl::json::parse( omemfile.str() );
    }
    catch ( nl::json::parse_error& e )
    {
        std::cout << "Error parsing JSON response: " << e.what() << "\n";
        return {};
    }

    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( jsonResponse, "Create token fails", status.code() ) << "\n";
        return {};
    }

    std::string tokenCreated;
    if ( jsonResponse.contains( "authToken" ) && jsonResponse["authToken"].contains( "token" ) )
    {
        tokenCreated = jsonResponse["authToken"]["token"].get<std::string>();
    }
    if (!validateToken(tokenCreated))
    {
        std::cout << Girder::errorMessage( jsonResponse, "Validation token fails", status.code() ) << "\n";
        return {};
    }
    return tokenCreated;
}

bool RemoteData::Girder::validateToken(const std::string& token) const
{
    try
    {
        std::cout << fmt::format("Validating token: {}", token) << "\n";
        std::string url = M_url + "/api/v1/user/me";
        cpr::Response res = cpr::Get(
            cpr::Url{url},
            cpr::Header{{"Girder-Token", token}},
            cpr::VerifySsl{false}
        );

        if (res.status_code == 200)
        {
            std::cout << "Token is valid.\n";
            return true;
        }
        else
        {
            std::cout << "Token validation failed: " << res.text << "\n";
            return false;
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "Exception during token validation: " << e.what() << "\n";
        return false;
    }
}
void RemoteData::Girder::removeToken( std::string const& token ) const
{
    std::string urlRemoveToken = M_url + "/api/v1/token/session";
    std::vector<std::string> headers{
        "Accept: application/json",
        "Girder-Token: " + token };

    std::ostringstream omemfile;
    StatusRequestHTTP status = requestHTTPCUSTOM( "DELETE", urlRemoveToken, headers, omemfile );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPCUSTOM (DELETE): " << status.msg() << "\n";
        return;
    }
    if ( status.code() != 200 )
    {
        // Parse JSON response
        nl::json jsonResponse;
        try
        {
            jsonResponse = nl::json::parse( omemfile.str() );
        }
        catch ( nl::json::parse_error& e )
        {
            std::cout << "Error parsing JSON response: " << e.what() << "\n";
            return;
        }

        std::cout << Girder::errorMessage( jsonResponse, "Remove token fails", status.code() ) << "\n";
        return;
    }
}

std::tuple<std::vector<std::shared_ptr<RemoteData::FolderInfo>>, std::vector<std::shared_ptr<RemoteData::ItemInfo>>, std::vector<std::shared_ptr<RemoteData::FileInfo>>>
RemoteData::Girder::contents() const
{
    auto res = std::make_tuple( std::vector<std::shared_ptr<FolderInfo>>(),
                                std::vector<std::shared_ptr<ItemInfo>>(),
                                std::vector<std::shared_ptr<FileInfo>>() );
    if ( M_worldComm->isMasterRank() )
    {
        // use token if given else create token if api key given
        std::string token = M_token;
        if ( M_token.empty() && !M_apiKey.empty() )
            token = this->createToken();

        for ( std::string const& folderId : M_folderIds )
        {
            auto resFolder = folderContentsImpl( folderId, token );
            if ( resFolder )
                std::get<0>( res ).push_back( resFolder );
        }
        for ( std::string const& itemId : M_itemIds )
        {
            auto resItem = itemInfoImpl( itemId, token );
            this->updateFilesImpl( resItem, token );
            if ( resItem )
                std::get<1>( res ).push_back( resItem );
        }
        for ( std::string const& fileId : M_fileIds )
        {
            auto resFile = fileInfoImpl( fileId, token );
            if ( resFile )
                std::get<2>( res ).push_back( resFile );
        }

        // delete token if created
        if ( M_token.empty() && !M_apiKey.empty() && !token.empty() )
            this->removeToken( token );
    }
    // if ( sync )

    return res;
}

std::shared_ptr<RemoteData::FileInfo>
RemoteData::Girder::fileInfoImpl( std::string const& fileId, std::string const& token ) const
{
    std::shared_ptr<FileInfo> res;

    std::string url = M_url + "/api/v1/file/" + fileId;

    std::vector<std::string> headers{
        "Accept: application/json" };
    if ( !token.empty() )
        headers.push_back( "Girder-Token: " + token );

    std::ostringstream omemfile;
    StatusRequestHTTP status = requestHTTPGET( url, headers, omemfile );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPGET: " << status.msg() << "\n";
        return res;
    }

    // Parse JSON response
    nl::json jsonResponse;
    try
    {
        jsonResponse = nl::json::parse( omemfile.str() );
    }
    catch ( nl::json::parse_error& e )
    {
        std::cout << "Error parsing JSON response: " << e.what() << "\n";
        return res;
    }

    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( jsonResponse, "File info fails", status.code() ) << "\n";
        return res;
    }

    std::string name = jsonResponse.value( "name", "" );
    std::string id = jsonResponse.value( "_id", "" );
    std::string mimeType = jsonResponse.value( "mimeType", "" );
    std::string sha512 = jsonResponse.value( "sha512", "" );
    size_type size = jsonResponse.value( "size", invalid_v<size_type> );

    res = std::make_shared<FileInfo>( name, id, size );
    res->setMimeType( mimeType );
    if ( !sha512.empty() )
        res->setChecksum( "sha512", sha512 );
    return res;
}
std::shared_ptr<RemoteData::FolderInfo>
RemoteData::Girder::folderInfoImpl( std::string const& folderId, std::string const& token ) const
{
    std::shared_ptr<FolderInfo> res;

    std::string url = M_url + "/api/v1/folder/" + folderId;

    std::vector<std::string> headers{
        "Accept: application/json" };
    if ( !token.empty() )
        headers.push_back( "Girder-Token: " + token );

    std::ostringstream omemfile;
    StatusRequestHTTP status = requestHTTPGET( url, headers, omemfile );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPGET: " << status.msg() << "\n";
        return res;
    }

    // Parse JSON response
    nl::json jsonResponse;
    try
    {
        jsonResponse = nl::json::parse( omemfile.str() );
    }
    catch ( nl::json::parse_error& e )
    {
        std::cout << "Error parsing JSON response: " << e.what() << "\n";
        return res;
    }

    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( jsonResponse, "Folder info fails", status.code() ) << "\n";
        return res;
    }

    std::string name = jsonResponse.value( "name", "" );
    std::string id = jsonResponse.value( "_id", "" );
    size_type size = jsonResponse.value( "size", invalid_v<size_type> );

    res = std::make_shared<FolderInfo>( name, id, size );
    return res;
}
std::shared_ptr<RemoteData::ItemInfo>
RemoteData::Girder::itemInfoImpl( std::string const& itemId, std::string const& token ) const
{
    std::shared_ptr<ItemInfo> res;

    std::string url = M_url + "/api/v1/item/" + itemId;

    std::vector<std::string> headers{
        "Accept: application/json" };
    if ( !token.empty() )
        headers.push_back( "Girder-Token: " + token );

    std::ostringstream omemfile;
    StatusRequestHTTP status = requestHTTPGET( url, headers, omemfile );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPGET: " << status.msg() << "\n";
        return res;
    }

    // Parse JSON response
    nl::json jsonResponse;
    try
    {
        jsonResponse = nl::json::parse( omemfile.str() );
    }
    catch ( nl::json::parse_error& e )
    {
        std::cout << "Error parsing JSON response: " << e.what() << "\n";
        return res;
    }

    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( jsonResponse, "Item info fails", status.code() ) << "\n";
        return res;
    }

    std::string name = jsonResponse.value( "name", "" );
    std::string id = jsonResponse.value( "_id", "" );
    size_type size = jsonResponse.value( "size", invalid_v<size_type> );

    res = std::make_shared<ItemInfo>( name, id, size );
    return res;
}

std::shared_ptr<RemoteData::FolderInfo>
RemoteData::Girder::folderContentsImpl( std::string const& folderId, std::string const& token ) const
{
    auto res = folderInfoImpl( folderId, token );
    if ( !res )
        return res;

    std::string urlFolder = M_url + "/api/v1/folder?parentType=folder";
    urlFolder += "&parentId=" + folderId;
    urlFolder += "&sort=lowerName&sortdir=1";

    std::vector<std::string> headers{
        "Accept: application/json" };
    if ( !token.empty() )
        headers.push_back( "Girder-Token: " + token );

    // Get folders inside the folder
    std::ostringstream omemfileFolder;
    StatusRequestHTTP status = requestHTTPGET( urlFolder, headers, omemfileFolder );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPGET: " << status.msg() << "\n";
        return res;
    }

    // Parse JSON response for folders
    nl::json jsonResponseFolders;
    try
    {
        jsonResponseFolders = nl::json::parse( omemfileFolder.str() );
    }
    catch ( nl::json::parse_error& e )
    {
        std::cout << "Error parsing JSON response: " << e.what() << "\n";
        return res;
    }

    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( jsonResponseFolders, "Folder contents (folder) fails", status.code() ) << "\n";
        return res;
    }

    for ( const auto& item : jsonResponseFolders )
    {
        std::string name = item.value( "name", "" );
        std::string id = item.value( "_id", "" );
        size_type size = item.value( "size", invalid_v<size_type> );
        res->addFolder( std::make_shared<FolderInfo>( name, id, size ) );
    }

    // Get items inside the folder
    std::string urlItem = M_url + "/api/v1/item";
    urlItem += "?folderId=" + folderId;
    urlItem += "&limit=0&sort=lowerName&sortdir=1";

    std::ostringstream omemfileItem;
    status = requestHTTPGET( urlItem, headers, omemfileItem );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPGET: " << status.msg() << "\n";
        return res;
    }

    // Parse JSON response for items
    nl::json jsonResponseItems;
    try
    {
        jsonResponseItems = nl::json::parse( omemfileItem.str() );
    }
    catch ( nl::json::parse_error& e )
    {
        std::cout << "Error parsing JSON response: " << e.what() << "\n";
        return res;
    }

    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( jsonResponseItems, "Folder contents (item) fails", status.code() ) << "\n";
        return res;
    }

    for ( const auto& item : jsonResponseItems )
    {
        std::string name = item.value( "name", "" );
        std::string id = item.value( "_id", "" );
        size_type size = item.value( "size", invalid_v<size_type> );
        auto itemInfo = std::make_shared<ItemInfo>( name, id, size );
        this->updateFilesImpl( itemInfo, token );
        res->addItem( itemInfo );
    }

    return res;
}

void RemoteData::Girder::updateFilesImpl( std::shared_ptr<RemoteData::ItemInfo> itemInfo, std::string const& token ) const
{
    if ( !itemInfo )
        return;
    std::string itemId = itemInfo->id();
    if ( itemId.empty() )
        return;

    std::string url = M_url + "/api/v1/item/" + itemId + "/files?limit=0&sort=name&sortdir=1";

    std::vector<std::string> headers{
        "Accept: application/json" };
    if ( !token.empty() )
        headers.push_back( "Girder-Token: " + token );

    std::ostringstream omemfile;
    StatusRequestHTTP status = requestHTTPGET( url, headers, omemfile );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPGET: " << status.msg() << "\n";
        return;
    }

    // Parse JSON response
    nl::json jsonResponse;
    try
    {
        jsonResponse = nl::json::parse( omemfile.str() );
    }
    catch ( nl::json::parse_error& e )
    {
        std::cout << "Error parsing JSON response: " << e.what() << "\n";
        return;
    }

    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( jsonResponse, "File in item fails", status.code() ) << "\n";
        return;
    }

    for ( const auto& ptFile : jsonResponse )
    {
        std::string name = ptFile.value( "name", "" );
        std::string id = ptFile.value( "_id", "" );
        std::string mimeType = ptFile.value( "mimeType", "" );
        std::string sha512 = ptFile.value( "sha512", "" );
        size_type size = ptFile.value( "size", invalid_v<size_type> );

        auto fileInfo = std::make_shared<FileInfo>( name, id, size );
        fileInfo->setMimeType( mimeType );
        if ( !sha512.empty() )
            fileInfo->setChecksum( "sha512", sha512 );
        itemInfo->add( fileInfo );
    }
}

nl::json
RemoteData::Girder::resourceLookup( const std::string& path, const std::string& token ) const
{
    // Construct the URL
    std::string url = fmt::format( "{}/api/v1/resource/lookup", M_url);

    // Set up headers
    cpr::Header headers;
    headers["Accept"] = "application/json";
    if ( !token.empty() )
    {
        headers["Girder-Token"] = token;
    }
    //std::cout << fmt::format( "Resource lookup: {}, path; {}", url, path ) << std::endl;
    // Make the GET request
    cpr::Response res = cpr::Get(
        cpr::Url{url},
        cpr::Parameters{{"path", path}}
        // If authentication is needed:
        , headers
    );

    if ( res.status_code != 200 )
    {
        LOG(INFO) << fmt::format( "Failed to lookup resource: HTTP {}", res.status_code );
        return {};
    }

    // Parse the JSON response
    nl::json jsonResponse = nl::json::parse( res.text );

    return jsonResponse;
}
bool
RemoteData::Girder::deleteResource(const nl::json& resource, const std::string& _token) const
{
    // Use token if given else create token if API key given
    std::string token = M_token;
    if (M_token.empty())
        token = this->createToken();

    // Ensure the resource contains required fields
    if (!resource.contains("_id") || !resource.contains("_modelType"))
    {
        std::cout << "Invalid resource: missing _id or _modelType\n";
        return false;
    }

    // Get the resource ID and resource type
    std::string resourceId = resource["_id"].get<std::string>();
    std::string resourceType = resource["_modelType"].get<std::string>();

    // Construct the URL for deleting the resource
    std::string url = fmt::format("{}/api/v1/{}/{}", M_url, resourceType, resourceId);

    // Set up headers
    cpr::Header headers = {
        {"Accept", "application/json"},
        {"Girder-Token", token}
    };

    // Send the DELETE request
    cpr::Response res = cpr::Delete(
        cpr::Url{url},
        headers,
        cpr::VerifySsl{false} // Set to true in production
    );

    // Check the HTTP status code
    if (res.status_code != 200)
    {
        // Parse the error response
        nl::json jsonResponse;
        try
        {
            jsonResponse = nl::json::parse(res.text);
        }
        catch (nl::json::parse_error& e)
        {
            std::cout << "Error parsing JSON response: " << e.what() << "\n";
            return false;
        }

        // Output the error message
        std::cout << Girder::errorMessage(jsonResponse, "Delete resource fails", res.status_code) << "\n";
        return false;
    }

    if (M_token.empty() && !token.empty())
    {
        // Delete the token if created
        this->removeToken(token);
    }
    // Successfully deleted the resource
    std::cout << "Resource deleted successfully: " << resourceId << " (" << resourceType << ")\n";
    return true;
}

} // namespace Feel