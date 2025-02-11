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

// Function to preprocess JSON-like input
inline std::string preprocessCustomFormat(const std::string& input)
{
    std::string output = input;

    // Add quotes around keys
    std::regex keyRegex(R"((\b[a-zA-Z_][a-zA-Z0-9_]*\b)\s*:)");
    output = std::regex_replace(output, keyRegex, R"("$1":)");

    // Add quotes around string values (excluding those in arrays or objects)
    std::regex valueRegex(R"(:\s*([a-zA-Z0-9_./\-]+)\s*[,\}])");
    output = std::regex_replace(output, valueRegex, R"(: "$1"$&)");

    // Fix arrays to keep them valid JSON
    std::regex arrayRegex(R"(":?\s*"\[.*?\]")");

    std::string temp_output;
    std::string::const_iterator start = output.begin();
    std::string::const_iterator end = output.end();

    std::regex_iterator<std::string::const_iterator> rit(start, end, arrayRegex);
    std::regex_iterator<std::string::const_iterator> rend;

    std::size_t last_pos = 0;

    for (; rit != rend; ++rit)
    {
        // Append text before the match
        temp_output.append(start + last_pos, rit->prefix().second);

        // Get the matched string
        std::string matchedStr = rit->str();

        // Remove the surrounding quotes
        size_t firstQuote = matchedStr.find('\"');
        size_t lastQuote = matchedStr.rfind('\"');

        if (firstQuote != std::string::npos && lastQuote != std::string::npos && firstQuote != lastQuote)
        {
            temp_output.append(matchedStr.substr(0, firstQuote)); // Before the first quote
            temp_output.append(matchedStr.substr(firstQuote + 1, lastQuote - firstQuote - 1)); // Inside quotes
            temp_output.append(matchedStr.substr(lastQuote + 1)); // After the last quote
        }
        else
        {
            // No quotes found, append the matched string as is
            temp_output.append(matchedStr);
        }

        // Update last position
        last_pos = rit->position() + rit->length();
    }

    // Append any remaining text after the last match
    temp_output.append(start + last_pos, end);

    output = temp_output;

    return output;
}
nl::json parseCustomFormat(const std::string& customInput)
{
    try
    {
        // Preprocess the input to valid JSON
        std::string jsonText = preprocessCustomFormat(customInput);

        // Parse using nlohmann::json
        return nl::json::parse(jsonText);
    }
    catch (const nl::json::exception& e)
    {
        std::cerr << "JSON parsing error: " << e.what() << std::endl;
        throw;
    }
}

StatusRequestHTTP requestHTTPGET( const std::string& url, const std::vector<std::string>& headers, std::ostream& ofile, int timeout, int max_retries, int backoff_delay  )
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
                                   std::ostream& ofile, int timeout, int max_retries, int backoff_delay )
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
                                   int timeout, int max_retries, int backoff_delay)
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
StatusRequestHTTP requestHTTPCUSTOM( const std::string& customRequest, const std::string& url, const std::vector<std::string>& headers, std::ostream& ofile, int timeout, int max_retries, int backoff_delay )
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

StatusRequestHTTP requestDownloadURL( const std::string& url, std::ostream& ofile, int timeout, int max_retries, int backoff_delay )
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
    //RemoteData::CKAN ckanTool( desc, *worldComm );
    //if ( ckanTool.isInit() )
    //{
    //    M_ckan.emplace( ckanTool );
    //    return;
    //}
}

bool RemoteData::canDownload() const
{
    if ( M_url )
        return true;
    else if ( M_github )
        return true;
    else if ( M_girder && M_girder->canDownload() )
        return true;
    else if ( M_ckan && M_ckan->canDownload() )
        return true;
    return false;
}

bool RemoteData::canUpload() const
{
    if ( M_girder && M_girder->canUpload() )
        return true;
    else if ( M_ckan && M_ckan->canUpload() )
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
    else if ( M_ckan )
        return M_ckan->download( dir );
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

// In RemoteData.hpp or RemoteData.cpp

nl::json RemoteData::createDataset(const std::string& datasetName) const
{
    if (M_ckan && M_ckan->isInit())
        return M_ckan->createDataset(datasetName);
    return nl::json(); // Return an empty JSON object if CKAN is not initialized
}

bool RemoteData::deleteDataset(const std::string& datasetId) const
{
    if (M_ckan && M_ckan->isInit())
        return M_ckan->deleteDataset(datasetId);
    return false; // Return false if CKAN is not initialized
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
