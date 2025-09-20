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
RemoteData::Github::Github( std::string const& desc, WorldComm& worldComm )
    : M_worldComm( worldComm.shared_from_this() )
{
    std::regex ex("([ ]*)github([ ]*):([ ]*)([{])([^]*)([}])");
    std::cmatch what;
    if( !regex_match(desc.c_str(), what, ex) )
        return;

    CHECK( what.size() == 7 ) << "invalid size";
    std::vector<std::string> keysvalues;
    std::string exprtosplit = std::string(what[5].first, what[5].second);


    auto resConvertion = convertDescToJson( exprtosplit );

    if ( !resConvertion.first )
    {
        if ( M_worldComm->isMasterRank() )
            std::cout << "Github desc has a syntax error: " << desc << "\n";
        return;
    }
    nl::json jsonObj = resConvertion.second;

    if ( jsonObj.contains( "owner" ) )
        M_owner = jsonObj["owner"].get<std::string>();
    if ( jsonObj.contains( "repo" ) )
        M_repo = jsonObj["repo"].get<std::string>();
    if ( jsonObj.contains( "branch" ) )
        M_branch = jsonObj["branch"].get<std::string>();
    if ( jsonObj.contains( "path" ) )
    {
        std::string path = jsonObj["path"].get<std::string>();
        if ( Feel::filename_is_dot( fs::path( path ).filename() ) )
            M_path = fs::path( path ).parent_path().string();
        else
            M_path = path;
    }
    if ( jsonObj.contains( "token" ) )
        M_token = jsonObj["token"].get<std::string>();
        // If token is empty, try looking for it in environment variable FEELPP_GITHUB_TOKEN
    if ( M_token.empty() )
    {
        char* env = std::getenv( "FEELPP_GITHUB_TOKEN" );
        if ( env )
        {
            M_token = env;
            std::cout << "Token loaded from environment variable" << std::endl;
        }
        else
        {
            std::cout << "No FEELPP_GITHUB_TOKEN environment variable found" << std::endl;
        }
    }
    else
    {
        std::cout << "Token loaded from configuration" << std::endl;
    }
    std::cout << "Token status: " << (M_token.empty() ? "empty" : "set") << std::endl;
    if ( M_owner.empty() )
        M_owner = "feelpp";
    if ( M_repo.empty() )
        M_repo = "feelpp";
}
bool RemoteData::Github::isInit() const
{
    return !M_owner.empty() && !M_repo.empty();
}
std::vector<std::string>
RemoteData::Github::download( std::string const& dir ) const
{
    std::vector<std::string> downloadFileOrFolder;
    if ( M_worldComm->isMasterRank() )
    {
        downloadFileOrFolder = this->downloadImpl( dir );
    }
    mpi::broadcast( M_worldComm->globalComm(), downloadFileOrFolder, M_worldComm->masterRank() );
    return downloadFileOrFolder;
}

std::vector<std::string>
RemoteData::Github::downloadImpl( std::string const& dir ) const
{
    std::vector<std::string> downloadFileOrFolder;

    // Determine progress level
    RemoteDataProgress::Level progressLevel = RemoteDataProgress::Level::NORMAL;
    if (Environment::vm().count("quiet")) {
        progressLevel = RemoteDataProgress::Level::QUIET;
    } else if (Environment::vm().count("debug")) {
        progressLevel = RemoteDataProgress::Level::DEBUG;
    } else if (Environment::vm().count("verbose") || Environment::vm().count("progress")) {
        progressLevel = RemoteDataProgress::Level::VERBOSE;
    }
    
    RemoteDataProgress progress(RemoteDataProgress::Operation::DOWNLOAD, progressLevel);

    std::string url = "https://api.github.com/repos/" + M_owner + "/" + M_repo + "/contents/" + M_path;
    if ( !M_branch.empty() )
        url += "?ref=" + M_branch;
    std::vector<std::string> headers;
    headers.push_back( "Accept: application/vnd.github.v3.json" );
    headers.push_back( "User-Agent: feelpp-agent" );
    // Temporarily disable token authentication to test if API works without it
    /*
    if ( !M_token.empty() )
    {
        headers.push_back( "Authorization: Bearer " + M_token );
        std::cout << "Using Bearer token for authentication" << std::endl;
    }
    else
    {
        std::cout << "No token available for authentication" << std::endl;
    }
    */
    if (progress.isDebug()) {
        progress.debug("Testing without authentication token");
    }

    // Progress and debug output
    std::string repo_info = M_owner + "/" + M_repo + (M_path.empty() ? "" : "/" + M_path);
    if (!M_branch.empty()) {
        repo_info += " (branch: " + M_branch + ")";
    }
    progress.startOperation("GitHub", repo_info);
    
    if (progress.isDebug()) {
        progress.debug("GitHub URL: " + url);
        progress.debug("Owner: " + M_owner + ", Repo: " + M_repo + ", Path: " + M_path + ", Branch: " + M_branch);
    }

    std::ostringstream omemfile;
    StatusRequestHTTP status = requestHTTPGET( url, headers, omemfile );
    if ( !status.success() )
    {
        progress.error("GitHub error in requestHTTPGET: " + status.msg());
        if (progress.isDebug()) {
            progress.debug("HTTP Status Code: " + std::to_string(status.code()));
        }
        return {};
    }

    // Parse the JSON response using nl::json
    nl::json jsonResponse;
    try
    {
        jsonResponse = nl::json::parse( omemfile.str() );
    }
    catch ( nl::json::parse_error& e )
    {
        progress.error("Error parsing JSON response: " + std::string(e.what()));
        return {};
    }

    if ( status.code() != 200 )
    {
        progress.error(Github::errorMessage( jsonResponse, "Get metadata fails", status.code() ));
        return {};
    }

    if ( jsonResponse.contains( "path" ) )
    {
        // It's a file
        std::string filename = jsonResponse["name"].get<std::string>();
        headers[0] = "Accept: application/vnd.github.v3.raw";
        std::string filepath = ( fs::path( dir ) / filename ).string();
        
        if (progress.isVerbose() || progress.isDebug()) {
            progress.startFile(filename, 0); // We don't know size from GitHub API without additional request
        }
        
        std::ofstream ofile( filepath, std::ios::out | std::ios::binary );
        status = requestHTTPGET( url, headers, ofile );
        ofile.close();
        if ( !status.success() )
        {
            progress.error("GitHub error in requestHTTPGET: " + status.msg());
            return {};
        }
        if ( status.code() != 200 )
        {
            progress.error(Github::errorMessage( nl::json{}, "Download file fails", status.code() ));
            return {};
        }

        if (progress.isNormal() || progress.isVerbose()) {
            // Get actual file size after download
            std::streamsize fileSize = 0;
            if (fs::exists(filepath)) {
                fileSize = fs::file_size(filepath);
            }
            progress.completeFile(filename, fileSize);
        }
        
        downloadFileOrFolder.push_back( filepath );
    }
    else
    {
        // It's a directory
        std::string subdir = ( M_path.empty() ) ? M_repo : fs::path( M_path ).filename().string();
        std::string newdir = ( fs::path( dir ) / subdir ).string();
        fs::create_directories( newdir );

        if (progress.isVerbose() || progress.isDebug()) {
            progress.debug("Downloading directory: " + subdir + " (" + std::to_string(jsonResponse.size()) + " items)");
        }

        auto resFolder = this->downloadFolderRecursively( jsonResponse, newdir, progress );
        if ( !std::get<0>( resFolder ) )
        {
            progress.error(std::get<1>( resFolder ));
            return {};
        }
        downloadFileOrFolder.push_back( newdir );
    }

    if (progress.isNormal() || progress.isVerbose()) {
        progress.completeOperation();
    }

    return downloadFileOrFolder;
}

std::tuple<bool, std::string>
RemoteData::Github::downloadFolderRecursively( nl::json const& jsonResponse, std::string const& dir, const RemoteDataProgress& progress ) const
{
    std::vector<std::string> headers;
    headers.push_back( "Accept: application/vnd.github.v3.json" );
    headers.push_back( "User-Agent: feelpp-agent" );
    if ( !M_token.empty() )
        headers.push_back( "Authorization: Bearer " + M_token );

    int itemCount = 0;
    int totalItems = jsonResponse.size();

    for ( const auto& item : jsonResponse )
    {
        itemCount++;
        std::string type = item["type"].get<std::string>();
        std::string name = item["name"].get<std::string>();
        std::string pathInUrl = item["path"].get<std::string>();
        std::string url = "https://api.github.com/repos/" + M_owner + "/" + M_repo + "/contents/" + pathInUrl;
        if ( !M_branch.empty() )
            url += "?ref=" + M_branch;

        if ( type == "file" )
        {
            // Download file
            headers[0] = "Accept: application/vnd.github.v3.raw";
            std::string filepath = ( fs::path( dir ) / name ).string();
            
            if (progress.isVerbose() || progress.isDebug()) {
                progress.startFile(name, 0, itemCount, totalItems); // We don't know size from GitHub API
            }
            
            std::ofstream ofile( filepath, std::ios::out | std::ios::binary );
            StatusRequestHTTP status = requestHTTPGET( url, headers, ofile );
            ofile.close();
            if ( !status.success() )
                return std::make_tuple( false, "GitHub error in requestHTTPGET : " + status.msg() );
            if ( status.code() != 200 )
                return std::make_tuple( false, Github::errorMessage( nl::json{}, "Download file fails", status.code() ) );
                
            if (progress.isNormal() || progress.isVerbose()) {
                // Get actual file size after download
                std::streamsize fileSize = 0;
                if (fs::exists(filepath)) {
                    fileSize = fs::file_size(filepath);
                }
                progress.completeFile(name, fileSize);
            }
        }
        else if ( type == "dir" )
        {
            // Get JSON for subdirectory
            headers[0] = "Accept: application/vnd.github.v3.json";
            std::ostringstream omemfile;
            StatusRequestHTTP status = requestHTTPGET( url, headers, omemfile );
            if ( !status.success() )
                return std::make_tuple( false, "GitHub error in requestHTTPGET : " + status.msg() );
            nl::json subdirJson = nl::json::parse( omemfile.str() );
            if ( status.code() != 200 )
                return std::make_tuple( false, Github::errorMessage( subdirJson, "Get metadata fails", status.code() ) );
            // Create subdirectory
            std::string newdir = ( fs::path( dir ) / name ).string();
            fs::create_directories( newdir );
            
            if (progress.isVerbose() || progress.isDebug()) {
                progress.debug("Processing directory: " + name + " (" + std::to_string(subdirJson.size()) + " items)");
            }
            
            // Recursive call
            auto resRecur = this->downloadFolderRecursively( subdirJson, newdir, progress );
            if ( !std::get<0>( resRecur ) )
                return std::make_tuple( false, std::get<1>( resRecur ) );
        }
    }
    return std::make_tuple( true, "" );
}
std::string
RemoteData::Github::errorMessage( nl::json const& jsonResponse, std::string const& defaultMsg, uint16_type statusCode )
{
    std::string errMsg = defaultMsg;
    if ( jsonResponse.contains( "message" ) )
        errMsg = jsonResponse["message"].get<std::string>();
    std::string docUrl;
    if ( jsonResponse.contains( "documentation_url" ) )
        docUrl = jsonResponse["documentation_url"].get<std::string>();
    std::string errMsgInfo;
    if ( statusCode != invalid_uint16_type_value )
        errMsgInfo += ( boost::format( "code %1%" ) % statusCode ).str();
    std::string res = "GitHub error";
    if ( !errMsgInfo.empty() )
        res += " (" + errMsgInfo + ")";
    if ( !errMsg.empty() )
        res += " : " + errMsg;
    if ( !docUrl.empty() )
        res += ( boost::format( " [ doc_url: %1% ]" ) % docUrl ).str();
    return res;
}

} // namespace Feel