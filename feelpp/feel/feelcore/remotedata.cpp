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

#include <boost/algorithm/string.hpp>
#include <cpr/cpr.h>
#include <feel/feelcore/remotedata.hpp>
#include <fstream>
#include <regex>

namespace Feel
{

// database taken from : https://code.msdn.microsoft.com/windowsapps/Get-Mimetype-from-a-cd7890af
const std::map<std::string, std::string> mimeTypes =
    {
        //{"***",    "application/octet-stream"},
        { ".csv", "text/csv" },
        { ".tsv", "text/tab-separated-values" },
        { ".tab", "text/tab-separated-values" },
        { ".html", "text/html" },
        { ".htm", "text/html" },
        { ".doc", "application/msword" },
        { ".docx", "application/vnd.openxmlformats-officedocument.wordprocessingml.document" },
        { ".ods", "application/x-vnd.oasis.opendocument.spreadsheet" },
        { ".odt", "application/vnd.oasis.opendocument.text" },
        { ".rtf", "application/rtf" },
        { ".sxw", "application/vnd.sun.xml.writer" },
        { ".txt", "text/plain" },
        { ".xls", "application/vnd.ms-excel" },
        { ".xlsx", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet" },
        { ".pdf", "application/pdf" },
        { ".ppt", "application/vnd.ms-powerpoint" },
        { ".pps", "application/vnd.ms-powerpoint" },
        { ".pptx", "application/vnd.openxmlformats-officedocument.presentationml.presentation" },
        { ".wmf", "image/x-wmf" },
        { ".atom", "application/atom+xml" },
        { ".xml", "application/xml" },
        { ".json", "application/json" },
        { ".js", "application/javascript" },
        { ".ogg", "application/ogg" },
        { ".ps", "application/postscript" },
        { ".woff", "application/x-woff" },
        { ".xhtml", "application/xhtml+xml" },
        { ".xht", "application/xhtml+xml" },
        { ".zip", "application/zip" },
        { ".gz", "application/x-gzip" },
        { ".rar", "application/rar" },
        { ".rm", "application/vnd.rn-realmedia" },
        { ".rmvb", "application/vnd.rn-realmedia-vbr" },
        { ".swf", "application/x-shockwave-flash" },
        { ".au", "audio/basic" },
        { ".snd", "audio/basic" },
        { ".mid", "audio/mid" },
        { ".rmi", "audio/mid" },
        { ".mp3", "audio/mpeg" },
        { ".aif", "audio/x-aiff" },
        { ".aifc", "audio/x-aiff" },
        { ".aiff", "audio/x-aiff" },
        { ".m3u", "audio/x-mpegurl" },
        { ".ra", "audio/vnd.rn-realaudio" },
        { ".ram", "audio/vnd.rn-realaudio" },
        { ".wav", "audio/x-wave" },
        { ".wma", "audio/x-ms-wma" },
        { ".m4a", "audio/x-m4a" },
        { ".bmp", "image/bmp" },
        { ".gif", "image/gif" },
        { ".jpe", "image/jpeg" },
        { ".jpeg", "image/jpeg" },
        { ".jpg", "image/jpeg" },
        { ".jfif", "image/jpeg" },
        { ".png", "image/png" },
        { ".svg", "image/svg+xml" },
        { ".tif", "image/tiff" },
        { ".tiff", "image/tiff" },
        { ".ico", "image/vnd.microsoft.icon" },
        { ".css", "text/css" },
        { ".bas", "text/plain" },
        { ".c", "text/plain" },
        { ".h", "text/plain" },
        { ".rtx", "text/richtext" },
        { ".mp2", "video/mpeg" },
        { ".mpa", "video/mpeg" },
        { ".mpe", "video/mpeg" },
        { ".mpeg", "video/mpeg" },
        { ".mpg", "video/mpeg" },
        { ".mpv2", "video/mpeg" },
        { ".mov", "video/quicktime" },
        { ".qt", "video/quicktime" },
        { ".lsf", "video/x-la-asf" },
        { ".lsx", "video/x-la-asf" },
        { ".asf", "video/x-ms-asf" },
        { ".asr", "video/x-ms-asf" },
        { ".asx", "video/x-ms-asf" },
        { ".avi", "video/x-msvideo" },
        { ".3gp", "video/3gpp" },
        { ".3gpp", "video/3gpp" },
        { ".3g2", "video/3gpp2" },
        { ".movie", "video/x-sgi-movie" },
        { ".mp4", "video/mp4" },
        { ".wmv", "video/x-ms-wmv" },
        { ".webm", "video/webm" },
        { ".m4v", "video/x-m4v" },
        { ".flv", "video/x-flv" } };

static size_t write_data( char /*void*/* ptr, size_t size, size_t nmemb, void* stream )
{
    // size_t written = fwrite(ptr, size, nmemb, (FILE *)stream);
    //((std::ofstream * ) stream)->write(ptr,nmemb);
    ( (std::ostream*)stream )->write( ptr, nmemb );
    return nmemb /*written*/;
}

static size_t read_data( char /*void*/* ptr, size_t size, size_t nmemb, void* stream )
{
    ( (std::istream*)stream )->read( ptr, nmemb );
    return nmemb /*written*/;
}

class StatusRequestHTTP : public std::tuple<bool, uint16_type, std::string>
{
    typedef std::tuple<bool, uint16_type, std::string> super_type;

  public:
    StatusRequestHTTP( bool s, std::string const& m = "" )
        : super_type( s, invalid_uint16_type_value, m )
    {
    }
    StatusRequestHTTP( bool s, uint16_type c, std::string const& m = "" )
        : super_type( s, c, m )
    {
    }
    StatusRequestHTTP( StatusRequestHTTP const& ) = default;
    StatusRequestHTTP( StatusRequestHTTP&& ) = default;
    StatusRequestHTTP& operator=( StatusRequestHTTP const& ) = default;

    bool success() const { return std::get<0>( *this ); }
    uint16_type code() const { return std::get<1>( *this ); }
    std::string const& msg() const { return std::get<2>( *this ); }
};

StatusRequestHTTP requestHTTPGET( const std::string& url, const std::vector<std::string>& headers, std::ostream& ofile, int timeout = 5000, int max_retries = 3, int backoff_delay = 1000 )
{
    // Convert headers to cpr::Header
    cpr::Header cpr_headers;
    for ( const auto& header : headers )
    {
        auto pos = header.find( ": " );
        if ( pos != std::string::npos )
        {
            cpr_headers[header.substr( 0, pos )] = header.substr( pos + 2 );
        }
    }

    int retries = 0;
    while ( retries <= max_retries )
    {
        // Perform the GET request with timeout
        auto response = cpr::Get( cpr::Url{ url }, cpr_headers, cpr::Timeout{ timeout } );

        // Check if the request was successful
        if ( response.status_code == 200 )
        {
            ofile << response.text;
            return StatusRequestHTTP( true, response.status_code, "" );
        }

        // Retry on certain failure conditions (e.g., timeout or server errors)
        if ( response.error.code == cpr::ErrorCode::OPERATION_TIMEDOUT || response.status_code >= 500 )
        {
            ++retries;
            if ( retries > max_retries )
            {
                return StatusRequestHTTP( false, response.status_code, fmt::format( "Max retries reached {} for url {}", response.error.message, url ) );
            }
            // Wait before retrying
            std::this_thread::sleep_for( std::chrono::milliseconds( backoff_delay ) );
        }
        else
        {
            // Non-retryable error
            return StatusRequestHTTP( false, response.status_code, response.error.message );
        }
    }

    return StatusRequestHTTP( false, 0, fmt::format( "Unknown error occurred after retries to get url {}", url ) );
}
StatusRequestHTTP requestHTTPPOST( const std::string& url, const std::vector<std::string>& headers,
                                   std::ostream& ofile, int timeout = 5000, int max_retries = 3, int backoff_delay = 1000 )
{
    // Convert headers to cpr::Header
    cpr::Header cpr_headers;
    for ( const auto& header : headers )
    {
        auto pos = header.find( ": " );
        if ( pos != std::string::npos )
        {
            cpr_headers[header.substr( 0, pos )] = header.substr( pos + 2 );
        }
    }

    int retries = 0;
    while ( retries <= max_retries )
    {
        // Perform the POST request with timeout
        auto response = cpr::Post( cpr::Url{ url },
                                   cpr_headers,
                                   // cpr::Body{body},
                                   cpr::Timeout{ timeout } );

        // Check if the request was successful
        if ( response.status_code == 200 )
        {
            ofile << response.text;
            return StatusRequestHTTP( true, response.status_code, "" );
        }

        // Retry on certain failure conditions (e.g., timeout or server errors)
        if ( response.error.code == cpr::ErrorCode::OPERATION_TIMEDOUT || response.status_code >= 500 )
        {
            ++retries;
            if ( retries > max_retries )
            {
                return StatusRequestHTTP( false, response.status_code,
                                          fmt::format( "Max retries reached. Error: {} URL: {}", response.error.message, url ) );
            }
            // Wait before retrying
            std::this_thread::sleep_for( std::chrono::milliseconds( backoff_delay ) );
        }
        else
        {
            // Non-retryable error
            return StatusRequestHTTP( false, response.status_code, response.error.message );
        }
    }

    return StatusRequestHTTP( false, 0, fmt::format( "Unknown error occurred after retries for URL: {}", url ) );
}

StatusRequestHTTP requestHTTPPOST( const std::string& url, const std::vector<std::string>& headers,
                                   std::istream& ifile, int fsize, std::ostream& ofile,
                                   int timeout = 5000, int max_retries = 3, int backoff_delay = 1000 )
{
    // Convert headers to cpr::Header
    cpr::Header cpr_headers;
    for ( const auto& header : headers )
    {
        auto pos = header.find( ": " );
        if ( pos != std::string::npos )
        {
            cpr_headers[header.substr( 0, pos )] = header.substr( pos + 2 );
        }
    }

    // Read the input stream into a string using std::ostringstream
    std::ostringstream oss;
    oss << ifile.rdbuf();
    std::string body = oss.str();

    int retries = 0;
    while ( retries <= max_retries )
    {
        // Perform the POST request with timeout and body
        auto response = cpr::Post( cpr::Url{ url },
                                   cpr_headers,
                                   cpr::Body{ body },
                                   cpr::Timeout{ timeout } );

        // Write the response body to the output stream if the request was successful
        if ( response.status_code == 200 )
        {
            ofile.write( response.text.c_str(), response.text.size() );
            return StatusRequestHTTP( true, response.status_code, "" );
        }

        // Retry on timeout or server errors (5xx)
        if ( response.error.code == cpr::ErrorCode::OPERATION_TIMEDOUT || response.status_code >= 500 )
        {
            ++retries;
            if ( retries > max_retries )
            {
                return StatusRequestHTTP( false, response.status_code,
                                          fmt::format( "Max retries reached. Error: {} URL: {}", response.error.message, url ) );
            }
            std::this_thread::sleep_for( std::chrono::milliseconds( backoff_delay ) );
        }
        else
        {
            // Non-retryable error
            return StatusRequestHTTP( false, response.status_code, response.error.message );
        }
    }

    return StatusRequestHTTP( false, 0, fmt::format( "Unknown error occurred after retries for URL: {}", url ) );
}
StatusRequestHTTP requestHTTPCUSTOM( const std::string& customRequest, const std::string& url, const std::vector<std::string>& headers, std::ostream& ofile, int timeout = 5000, int max_retries = 3, int backoff_delay = 1000 )
{
    // Convert headers to cpr::Header
    cpr::Header cpr_headers;
    for ( const auto& header : headers )
    {
        auto pos = header.find( ": " );
        if ( pos != std::string::npos )
        {
            cpr_headers[header.substr( 0, pos )] = header.substr( pos + 2 );
        }
    }

    int retries = 0;
    while ( retries <= max_retries )
    {
        cpr::Response response;

        try
        {
            if ( customRequest == "PUT" )
            {
                response = cpr::Put( cpr::Url{ url }, cpr_headers, cpr::Timeout{ timeout } );
            }
            else if ( customRequest == "DELETE" )
            {
                response = cpr::Delete( cpr::Url{ url }, cpr_headers, cpr::Timeout{ timeout } );
            }
            else if ( customRequest == "PATCH" )
            {
                response = cpr::Patch( cpr::Url{ url }, cpr_headers, cpr::Timeout{ timeout } );
            }
            else if ( customRequest == "OPTIONS" )
            {
                response = cpr::Options( cpr::Url{ url }, cpr_headers, cpr::Timeout{ timeout } );
            }
            else
            {
                return StatusRequestHTTP( false, 400, "Unsupported HTTP method: " + customRequest );
            }

            // Check for success
            if ( response.status_code >= 200 && response.status_code < 300 )
            {
                ofile << response.text;
                return StatusRequestHTTP( true, response.status_code, "" );
            }

            // Retry logic for server errors or timeout
            if ( response.error.code == cpr::ErrorCode::OPERATION_TIMEDOUT || response.status_code >= 500 )
            {
                ++retries;
                if ( retries > max_retries )
                {
                    return StatusRequestHTTP( false, response.status_code, "Max retries reached for " + url );
                }
                std::this_thread::sleep_for( std::chrono::milliseconds( backoff_delay ) );
            }
            else
            {
                return StatusRequestHTTP( false, response.status_code, response.error.message );
            }
        }
        catch ( const std::exception& ex )
        {
            return StatusRequestHTTP( false, 0, std::string( "Exception: " ) + ex.what() );
        }
    }

    return StatusRequestHTTP( false, 0, "Unknown error occurred." );
}

StatusRequestHTTP requestDownloadURL( const std::string& url, std::ostream& ofile, int timeout = 5000, int max_retries = 3, int backoff_delay = 1000 )
{
    int retries = 0;
    while ( retries <= max_retries )
    {
        // Perform the GET request with timeout
        auto response = cpr::Get( cpr::Url{ url }, cpr::Timeout{ timeout } );

        // Check if the request was successful
        if ( response.status_code == 200 )
        {
            ofile << response.text; // Write the response body to the output stream
            return StatusRequestHTTP( true, response.status_code, "" );
        }

        // Retry on certain failure conditions (e.g., timeout or server errors)
        if ( response.error.code == cpr::ErrorCode::OPERATION_TIMEDOUT || response.status_code >= 500 )
        {
            ++retries;
            if ( retries > max_retries )
            {
                return StatusRequestHTTP( false, response.status_code, fmt::format( "Max retries reached. Error: {}", response.error.message ) );
            }
            // Wait before retrying
            std::this_thread::sleep_for( std::chrono::milliseconds( backoff_delay ) );
        }
        else
        {
            // Non-retryable error
            return StatusRequestHTTP( false, response.status_code, response.error.message );
        }
    }

    return StatusRequestHTTP( false, 0, "Unknown error occurred after retries." );
}

RemoteData::RemoteData( std::string const& desc, worldcomm_ptr_t const& worldComm )
    : M_worldComm( worldComm )
{
    RemoteData::URL urlTool( desc, *worldComm );
    if ( urlTool.isValid() )
    {
        M_url.emplace( urlTool );
        return;
    }
    RemoteData::Github githubTool( desc, *worldComm );
    if ( githubTool.isInit() )
    {
        M_github.emplace( githubTool );
        return;
    }
    RemoteData::Girder girderTool( desc, *worldComm );
    if ( girderTool.isInit() )
    {
        M_girder.emplace( girderTool );
        return;
    }
}

bool RemoteData::canDownload() const
{
    if ( M_url )
        return true;
    else if ( M_github )
        return true;
    else if ( M_girder && M_girder->canDownload() )
        return true;
    return false;
}

bool RemoteData::canUpload() const
{
    if ( M_girder && M_girder->canUpload() )
        return true;
    return false;
}

std::vector<std::string>
RemoteData::download( std::string const& dir, std::string const& filename ) const
{
    std::vector<std::string> downloadedData;
    if ( M_url )
        downloadedData.push_back( M_url->download( dir, filename ) );
    else if ( M_github )
        return M_github->download( dir );
    else if ( M_girder )
        return M_girder->download( dir );
    return downloadedData;
}

nl::json
RemoteData::resourceLookup( std::string const& path, std::string const& token ) const
{
    if ( M_girder )
        return M_girder->resourceLookup( path, token );
    return {};
}
bool
RemoteData::deleteResource( nl::json const& resourceId, std::string const& token ) const
{
    if ( M_girder )
        return M_girder->deleteResource( resourceId, token );
    return false;
}
std::vector<std::string>
RemoteData::upload( std::string const& dataPath, std::string const& parentId, bool sync ) const
{
    if ( M_girder && M_girder->canUpload() )
        return M_girder->upload( dataPath, parentId, sync );
    return {};
}

std::vector<std::vector<std::string>>
RemoteData::upload( std::vector<std::pair<std::string, std::string>> const& dataToUpload, bool sync ) const
{
    if ( M_girder && M_girder->canUpload() )
        return M_girder->upload( dataToUpload, sync );
    return {};
}

void RemoteData::replaceFile( std::string const& filePath, std::string const& fileId ) const
{
    if ( M_girder && M_girder->canUpload() )
        return M_girder->replaceFile( filePath, fileId );
}

void RemoteData::replaceFile( std::vector<std::pair<std::string, std::string>> const& filesToReplace ) const
{
    if ( M_girder && M_girder->canUpload() )
        return M_girder->replaceFile( filesToReplace );
}

std::vector<std::pair<std::string, std::string>>
RemoteData::createFolder( std::string const& folderPath, std::string const& parentId, bool sync ) const
{
    if ( M_girder && M_girder->canUpload() )
        return M_girder->createFolder( folderPath, parentId, sync );
    return {};
}

std::vector<std::pair<std::string, std::string>>
RemoteData::createItem( std::string const& itemPath, std::string const& parentId, bool sync ) const
{
    if ( M_girder && M_girder->canUpload() )
        return M_girder->createItem( itemPath, parentId, sync );
    return {};
}

RemoteData::ContentsInfo
RemoteData::contents() const
{
    if ( M_girder && M_girder->isInit() )
        return ContentsInfo{ M_girder->contents() };
    return ContentsInfo{ std::make_tuple( std::vector<std::shared_ptr<FolderInfo>>(), std::vector<std::shared_ptr<ItemInfo>>(), std::vector<std::shared_ptr<FileInfo>>() ) };
}

RemoteData::URL::URL( std::string const& url, WorldComm& worldComm )
    : M_worldComm( worldComm.shared_from_this() )
{
    std::regex ex( "(http|https)://([^/ :]+):?([^/ ]*)(/?[^ #?]*)\\x3f?([^ #]*)#?([^ ]*)" );
    std::cmatch what;
    if ( regex_match( url.c_str(), what, ex ) )
    {
        M_url = url;
        M_protocol = std::string( what[1].first, what[1].second );
        M_domain = std::string( what[2].first, what[2].second );
        M_port = std::string( what[3].first, what[3].second );
        M_path = std::string( what[4].first, what[4].second );
        M_query = std::string( what[5].first, what[5].second );
#if 0
        std::cout << "[" << url << "]" << std::endl;
        std::cout << "protocol: " << M_protocol << std::endl;
        std::cout << "domain: " << M_domain << std::endl;
        std::cout << "port:" << M_port << std::endl;
        std::cout << "path:" << M_path << std::endl;
        std::cout << "query:" << M_query << std::endl;
#endif
    }
}

bool RemoteData::URL::isValid() const
{
    return !M_url.empty();
}

std::string
RemoteData::URL::download( std::string const& _dir, std::string const& _filename ) const
{
    CHECK( this->isValid() ) << "url is not valid";

    std::string filename = _filename;
    if ( filename.empty() )
    {
        fs::path p = fs::path( M_path );
        if ( p.has_filename() && !Feel::filename_is_dot( p ) && !Feel::filename_is_dot_dot( p ) )
            filename = p.filename().string();
        else
            filename = "download";
    }

    fs::path dir = fs::path( _dir );
    bool addSubDirHashURL = false;
    if ( addSubDirHashURL )
    {
        size_t hashURL = std::hash<std::string>()( M_url );
        std::string subdir = ( boost::format( "%1%" ) % hashURL ).str();
        dir /= fs::path( subdir );
    }
    std::string thefilename = ( dir / filename ).string();

    std::string url = M_url;

    if ( M_worldComm->isMasterRank() )
    {
        if ( !fs::exists( dir ) )
            fs::create_directories( dir );

        std::ofstream ofile( thefilename, std::ios::out | std::ios::binary );
        /* open the file */
        if ( ofile )
        {
            StatusRequestHTTP status = requestDownloadURL( url, ofile );
            if ( !status.success() )
                std::cout << "Download error : " << status.msg() << "\n";
            if ( status.code() != 200 ) // is it really true for all type http,ftp,... ???
                std::cout << "Download error : returned code " << status.code() << "\n";

            ofile.close();
        }
        else
            std::cout << "Download error : failure when create file  " << thefilename << "\n";
    }

    M_worldComm->barrier();
    return thefilename;
}

std::pair<bool, nl::json>
convertDescToJson( std::string const& desc )
{
    // split the desc with respect to special characters
    std::vector<std::pair<bool,std::string>> descSplitted;
    std::string currentSplit;
    std::vector<char> splitChars = { ':',',','{','}','[',']' };
    char lastSplitChar = '0';
    for ( std::string::const_iterator it=desc.begin(); it!=desc.end(); ++it)
    {
        char c = *it;
        auto itFindSplitChar = std::find(splitChars.begin(),splitChars.end(), c );
        if ( ( itFindSplitChar != splitChars.end() ) && ( c != ':' || lastSplitChar != ':' ) )
        {
            boost::trim( currentSplit );
            if ( !currentSplit.empty() )
                descSplitted.push_back( std::make_pair( false,currentSplit ) );
            descSplitted.push_back( std::make_pair( true,std::string(1,c) ) );
            lastSplitChar = c;
            currentSplit.clear();
        }
        else
            currentSplit += c;
    }
    boost::trim( currentSplit );
    if ( !currentSplit.empty() )
    {
        descSplitted.push_back( std::make_pair( false, currentSplit ) );
    }

    // create new string convertible to property tree (by adding double quote)
    std::string newDesc = "{";
    for (int k=0;k<descSplitted.size();++k)
    {
        std::string expr = descSplitted[k].second;
        if ( descSplitted[k].first )
        {
            newDesc += expr;
            continue;
        }
        newDesc += "\"" + expr + "\"";
    }
    newDesc += "}";

    // Parse the description string into a JSON object
    try
    {
        nl::json jsonObj = nl::json::parse( newDesc );
        return std::make_pair( true, jsonObj );
    }
    catch ( nl::json::parse_error& e )
    {
        std::cout << "Error parsing JSON description: " << e.what() << "\n";
        return std::make_pair( false, nl::json{} );
    }
}

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
        if ( env != NULL && env[0] != '\0' )
        {
            M_token = env;
        }
    }
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

    std::string url = "https://api.github.com/repos/" + M_owner + "/" + M_repo + "/contents/" + M_path;
    if ( !M_branch.empty() )
        url += "?ref=" + M_branch;
    std::vector<std::string> headers;
    headers.push_back( "Accept: application/vnd.github.v3.json" );
    headers.push_back( "User-Agent: feelpp-agent" );
    if ( !M_token.empty() )
        headers.push_back( "Authorization: token " + M_token );

    std::ostringstream omemfile;
    StatusRequestHTTP status = requestHTTPGET( url, headers, omemfile );
    if ( !status.success() )
    {
        std::cout << "Github error in requestHTTPGET: " << status.msg() << "\n";
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
        std::cout << "Error parsing JSON response: " << e.what() << "\n";
        return {};
    }

    if ( status.code() != 200 )
    {
        std::cout << Github::errorMessage( jsonResponse, "Get metadata fails", status.code() ) << "\n";
        return {};
    }

    if ( jsonResponse.contains( "path" ) )
    {
        // It's a file
        std::string filename = jsonResponse["name"].get<std::string>();
        headers[0] = "Accept: application/vnd.github.v3.raw";
        std::string filepath = ( fs::path( dir ) / filename ).string();
        std::ofstream ofile( filepath, std::ios::out | std::ios::binary );
        status = requestHTTPGET( url, headers, ofile );
        ofile.close();
        if ( !status.success() )
        {
            std::cout << "Github error in requestHTTPGET: " << status.msg() << "\n";
            return {};
        }
        if ( status.code() != 200 )
        {
            std::cout << Github::errorMessage( nl::json{}, "Download file fails", status.code() ) << "\n";
            return {};
        }

        downloadFileOrFolder.push_back( filepath );
    }
    else
    {
        // It's a directory
        std::string subdir = ( M_path.empty() ) ? M_repo : fs::path( M_path ).filename().string();
        std::string newdir = ( fs::path( dir ) / subdir ).string();
        fs::create_directories( newdir );

        auto resFolder = this->downloadFolderRecursively( jsonResponse, newdir );
        if ( !std::get<0>( resFolder ) )
        {
            std::cout << std::get<1>( resFolder ) << "\n";
            return {};
        }
        downloadFileOrFolder.push_back( newdir );
    }

    return downloadFileOrFolder;
}

std::tuple<bool, std::string>
RemoteData::Github::downloadFolderRecursively( nl::json const& jsonResponse, std::string const& dir ) const
{
    std::vector<std::string> headers;
    headers.push_back( "Accept: application/vnd.github.v3.json" );
    headers.push_back( "User-Agent: feelpp-agent" );
    if ( !M_token.empty() )
        headers.push_back( "Authorization: token " + M_token );

    for ( const auto& item : jsonResponse )
    {
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
            std::ofstream ofile( filepath, std::ios::out | std::ios::binary );
            StatusRequestHTTP status = requestHTTPGET( url, headers, ofile );
            ofile.close();
            if ( !status.success() )
                return std::make_tuple( false, "GitHub error in requestHTTPGET : " + status.msg() );
            if ( status.code() != 200 )
                return std::make_tuple( false, Github::errorMessage( nl::json{}, "Download file fails", status.code() ) );
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
            // Recursive call
            auto resRecur = this->downloadFolderRecursively( subdirJson, newdir );
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
                // It's a parentId, need to determine the resource type
                // For backward compatibility, assume it's a folder
                resourceType = "folder";
                resourceId = destination;
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
    return tokenCreated;
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
std::ostringstream
RemoteData::FolderInfo::print( size_t nTab ) const
{
    std::ostringstream ostr;
    std::string tab;
    for ( int k = 0; k < nTab; ++k )
        tab += "|   ";
    ostr << tab << "|-- " << this->name() << "  [type:folder";
    if ( !this->id().empty() )
        ostr << ", id:" << this->id();
    if ( this->size() != invalid_v<size_type> )
        ostr << ", size:" << this->size();
    ostr << "]\n";

    for ( auto const& folderInfo : M_folders )
        ostr << folderInfo->print( nTab + 1 ).str();
    for ( auto const& itemInfo : M_items )
        ostr << itemInfo->print( nTab + 1 ).str();
    for ( auto const& fileInfo : M_files )
        ostr << fileInfo->print( nTab + 1 ).str();
    return ostr;
}

std::ostringstream
RemoteData::ItemInfo::print( size_t nTab ) const
{
    std::ostringstream ostr;
    std::string tab;
    for ( int k = 0; k < nTab; ++k )
        tab += "|   ";
    ostr << tab << "|-- " << this->name() << "  [type:item";
    if ( !this->id().empty() )
        ostr << ", id:" << this->id();
    if ( this->size() != invalid_v<size_type> )
        ostr << ", size:" << this->size();
    ostr << "]\n";

    for ( auto const& fileInfo : M_files )
        ostr << fileInfo->print( nTab + 1, true ).str();
    return ostr;
}
std::ostringstream
RemoteData::FileInfo::print( size_t nTab, bool isFileInItem ) const
{
    std::ostringstream ostr;
    std::string tab;
    for ( int k = 0; k < nTab; ++k )
        tab += "|   ";
    ostr << tab;
    if ( !isFileInItem )
        ostr << "|-- ";
    else
        ostr << "* ";
    ostr << this->name() << "  [type:file";
    if ( !this->id().empty() )
        ostr << ", id:" << this->id();
    if ( this->size() != invalid_v<size_type> )
        ostr << ", size:" << this->size();
    ostr << "]\n";
    return ostr;
}

} // namespace Feel
