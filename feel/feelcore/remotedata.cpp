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

#include <feel/feelcore/remotedata.hpp>
#include <fstream>
#include <regex>
#include <boost/algorithm/string.hpp>

extern "C" {

#include <curl/curl.h>

}
namespace Feel
{

// database taken from : https://code.msdn.microsoft.com/windowsapps/Get-Mimetype-from-a-cd7890af
const std::map<std::string,std::string> mimeTypes =
{
    //{"***",    "application/octet-stream"},
    {".csv",    "text/csv"},
    {".tsv",    "text/tab-separated-values"},
    {".tab",    "text/tab-separated-values"},
    {".html",    "text/html"},
    {".htm",    "text/html"},
    {".doc",    "application/msword"},
    {".docx",    "application/vnd.openxmlformats-officedocument.wordprocessingml.document"},
    {".ods",    "application/x-vnd.oasis.opendocument.spreadsheet"},
    {".odt",    "application/vnd.oasis.opendocument.text"},
    {".rtf",    "application/rtf"},
    {".sxw",    "application/vnd.sun.xml.writer"},
    {".txt",    "text/plain"},
    {".xls",    "application/vnd.ms-excel"},
    {".xlsx",    "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"},
    {".pdf",    "application/pdf"},
    {".ppt",    "application/vnd.ms-powerpoint"},
    {".pps",    "application/vnd.ms-powerpoint"},
    {".pptx",    "application/vnd.openxmlformats-officedocument.presentationml.presentation"},
    {".wmf",    "image/x-wmf"},
    {".atom",    "application/atom+xml"},
    {".xml",    "application/xml"},
    {".json",    "application/json"},
    {".js",    "application/javascript"},
    {".ogg",    "application/ogg"},
    {".ps",    "application/postscript"},
    {".woff",    "application/x-woff"},
    {".xhtml","application/xhtml+xml"},
    {".xht",    "application/xhtml+xml"},
    {".zip",    "application/zip"},
    {".gz",    "application/x-gzip"},
    {".rar",    "application/rar"},
    {".rm",    "application/vnd.rn-realmedia"},
    {".rmvb",    "application/vnd.rn-realmedia-vbr"},
    {".swf",    "application/x-shockwave-flash"},
    {".au",        "audio/basic"},
    {".snd",    "audio/basic"},
    {".mid",    "audio/mid"},
    {".rmi",        "audio/mid"},
    {".mp3",    "audio/mpeg"},
    {".aif",    "audio/x-aiff"},
    {".aifc",    "audio/x-aiff"},
    {".aiff",    "audio/x-aiff"},
    {".m3u",    "audio/x-mpegurl"},
    {".ra",    "audio/vnd.rn-realaudio"},
    {".ram",    "audio/vnd.rn-realaudio"},
    {".wav",    "audio/x-wave"},
    {".wma",    "audio/x-ms-wma"},
    {".m4a",    "audio/x-m4a"},
    {".bmp",    "image/bmp"},
    {".gif",    "image/gif"},
    {".jpe",    "image/jpeg"},
    {".jpeg",    "image/jpeg"},
    {".jpg",    "image/jpeg"},
    {".jfif",    "image/jpeg"},
    {".png",    "image/png"},
    {".svg",    "image/svg+xml"},
    {".tif",    "image/tiff"},
    {".tiff",    "image/tiff"},
    {".ico",    "image/vnd.microsoft.icon"},
    {".css",    "text/css"},
    {".bas",    "text/plain"},
    {".c",        "text/plain"},
    {".h",        "text/plain"},
    {".rtx",    "text/richtext"},
    {".mp2",    "video/mpeg"},
    {".mpa",    "video/mpeg"},
    {".mpe",    "video/mpeg"},
    {".mpeg",    "video/mpeg"},
    {".mpg",    "video/mpeg"},
    {".mpv2",    "video/mpeg"},
    {".mov",    "video/quicktime"},
    {".qt",    "video/quicktime"},
    {".lsf",    "video/x-la-asf"},
    {".lsx",    "video/x-la-asf"},
    {".asf",    "video/x-ms-asf"},
    {".asr",    "video/x-ms-asf"},
    {".asx",    "video/x-ms-asf"},
    {".avi",    "video/x-msvideo"},
    {".3gp",    "video/3gpp"},
    {".3gpp",    "video/3gpp"},
    {".3g2",    "video/3gpp2"},
    {".movie","video/x-sgi-movie"},
    {".mp4",    "video/mp4"},
    {".wmv",    "video/x-ms-wmv"},
    {".webm","video/webm"},
    {".m4v",    "video/x-m4v"},
    {".flv",    "video/x-flv"}
};

static size_t write_data(char/*void*/ *ptr, size_t size, size_t nmemb, void *stream)
{
    //size_t written = fwrite(ptr, size, nmemb, (FILE *)stream);
    //((std::ofstream * ) stream)->write(ptr,nmemb);
    ((std::ostream * ) stream)->write(ptr,nmemb);
    return nmemb/*written*/;
}

static size_t read_data(char/*void*/ *ptr, size_t size, size_t nmemb, void *stream) {
    ((std::istream * ) stream)->read(ptr,nmemb);
    return nmemb/*written*/;

}

class StatusRequestHTTP : public std::tuple<bool,uint16_type,std::string>
{
    typedef std::tuple<bool,uint16_type,std::string> super_type;
public :
    StatusRequestHTTP( bool s, std::string const& m = "" )
        :
        super_type( s, invalid_uint16_type_value, m )
        {}
    StatusRequestHTTP( bool s, uint16_type c, std::string const& m = "" )
        :
        super_type( s, c, m )
        {}
    StatusRequestHTTP( StatusRequestHTTP const& ) = default;
    StatusRequestHTTP( StatusRequestHTTP && ) = default;
    StatusRequestHTTP& operator=( StatusRequestHTTP const& ) = default;

    bool success() const { return std::get<0>( *this ); }
    uint16_type code() const { return std::get<1>( *this ); }
    std::string const& msg() const { return std::get<2>( *this ); }
};

StatusRequestHTTP requestHTTPGET( std::string const& url, std::vector<std::string> const& headers, std::ostream & ofile )
{
#if defined(FEELPP_HAS_LIBCURL)
    CURLcode res;
    res = curl_global_init(CURL_GLOBAL_ALL);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );

    // init
    CURL *curl_handle = curl_easy_init();
    if ( !curl_handle ) return StatusRequestHTTP( false, "fail to run curl_easy_init" );

    // url
    res = curl_easy_setopt(curl_handle, CURLOPT_URL, url.c_str() );
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );
    // request GET type
    res = curl_easy_setopt(curl_handle, CURLOPT_HTTPGET, 1L);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );
    // headers
    struct curl_slist *list = NULL;
    for ( std::string const& header : headers )
        list = curl_slist_append(list, header.c_str() );
    res = curl_easy_setopt(curl_handle, CURLOPT_HTTPHEADER, list);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );

    /* send all data to this function  */
    res = curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, write_data);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );

    /* write the page body to this file handle */
    res = curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, &ofile);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );

    /* get it! */
    res = curl_easy_perform(curl_handle);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );

    // get http code status
    long http_code = 0;
    res = curl_easy_getinfo(curl_handle, CURLINFO_RESPONSE_CODE, &http_code);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );

    /* cleanup curl stuff */
    curl_easy_cleanup(curl_handle);
    curl_global_cleanup();

    return StatusRequestHTTP( true, http_code, "" );
#else
    CHECK( false ) << "LIBCURL is not detected";
    return StatusRequestHTTP( false );
#endif
}

StatusRequestHTTP requestHTTPPOST( std::string const& url, std::vector<std::string> const& headers, std::ostream & ofile )
{
#if defined(FEELPP_HAS_LIBCURL)
    CURLcode res;
    res = curl_global_init(CURL_GLOBAL_ALL);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );

    // init
    CURL *curl_handle = curl_easy_init();
    if ( !curl_handle ) return StatusRequestHTTP( false, "fail to run curl_easy_init" );
    //url
    res = curl_easy_setopt(curl_handle, CURLOPT_URL, url.c_str() );
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );
    // request POST type
    res = curl_easy_setopt(curl_handle, CURLOPT_POST, 1L);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );
    // headers
    struct curl_slist *list = NULL;
    for ( std::string const& header : headers )
        list = curl_slist_append(list, header.c_str() );
    res = curl_easy_setopt(curl_handle, CURLOPT_HTTPHEADER, list);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );

    // post an empty file
    long fsize = 0;
    char * postthis = nullptr;
    res = curl_easy_setopt(curl_handle, CURLOPT_POSTFIELDS, postthis);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );
    res = curl_easy_setopt(curl_handle, CURLOPT_POSTFIELDSIZE, fsize);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );

    /* send all data to this function  */
    res = curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, write_data);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );
    /* write the page body to this file handle */
    res = curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, &ofile);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );
    /* get it! */
    res = curl_easy_perform(curl_handle);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );

    // get http code status
    long http_code = 0;
    res = curl_easy_getinfo(curl_handle, CURLINFO_RESPONSE_CODE, &http_code);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );

    /* cleanup curl stuff */
    curl_easy_cleanup(curl_handle);
    curl_global_cleanup();

    return StatusRequestHTTP( true, http_code, "" );
#else
    CHECK( false ) << "LIBCURL is not detected";
    return StatusRequestHTTP( false );
#endif
}

void requestHTTPPOST( std::string const& url, std::vector<std::string> const& headers, std::istream & ifile, long fsize, std::ostream & ofile )
{
#if defined(FEELPP_HAS_LIBCURL)
    curl_global_init(CURL_GLOBAL_ALL);
    CURL *curl_handle;
    curl_handle = curl_easy_init();

    curl_easy_setopt(curl_handle, CURLOPT_URL, url.c_str() );

#if 1
    char * postthis = new char [fsize];
    ifile.read( postthis,fsize );
    curl_easy_setopt(curl_handle, CURLOPT_POSTFIELDS, postthis);
#else
    // very slow (need to understand)
    curl_easy_setopt(curl_handle, CURLOPT_POST, 1L);
    curl_easy_setopt(curl_handle, CURLOPT_READDATA, &ifile);
    curl_easy_setopt(curl_handle, CURLOPT_READFUNCTION, read_data);
#endif
    curl_easy_setopt(curl_handle, CURLOPT_POSTFIELDSIZE, fsize);

    struct curl_slist *list = NULL;
    for ( std::string const& header : headers )
        list = curl_slist_append(list, header.c_str() );
    curl_easy_setopt(curl_handle, CURLOPT_HTTPHEADER, list);

    /* send all data to this function  */
    curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, write_data);
    /* write the page body to this file handle */
    curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, &ofile);

    /* get it! */
    curl_easy_perform(curl_handle);

    /* cleanup curl stuff */
    curl_easy_cleanup(curl_handle);

    curl_global_cleanup();

    delete [] postthis;
#else
    CHECK( false ) << "LIBCURL is not detected";
#endif
}

void requestHTTPDELETE( std::string const& url, std::vector<std::string> const& headers, std::ostream & ofile )
{
#if defined(FEELPP_HAS_LIBCURL)
    curl_global_init(CURL_GLOBAL_ALL);
    CURL *curl_handle;
    curl_handle = curl_easy_init();

    curl_easy_setopt(curl_handle, CURLOPT_URL, url.c_str() );

    curl_easy_setopt(curl_handle, CURLOPT_CUSTOMREQUEST, "DELETE");

    struct curl_slist *list = NULL;
    for ( std::string const& header : headers )
        list = curl_slist_append(list, header.c_str() );
    curl_easy_setopt(curl_handle, CURLOPT_HTTPHEADER, list);

    /* send all data to this function  */
    curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, write_data);
    /* write the page body to this file handle */
    curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, &ofile);

    /* get it! */
    curl_easy_perform(curl_handle);

    /* cleanup curl stuff */
    curl_easy_cleanup(curl_handle);

    curl_global_cleanup();
#else
    CHECK( false ) << "LIBCURL is not detected";
#endif
}



RemoteData::RemoteData( std::string const& desc, WorldComm const& worldComm )
{
    RemoteData::URL urlTool( desc,worldComm );
    if ( urlTool.isValid() )
    {
        M_url.emplace( urlTool );
        return;
    }
    RemoteData::Github githubTool( desc,worldComm );
    if ( githubTool.isInit() )
    {
        M_github.emplace( githubTool );
        return;
    }
    RemoteData::Girder girderTool( desc,worldComm );
    if ( girderTool.isInit() )
    {
        M_girder.emplace( girderTool );
        return;
    }
}

bool
RemoteData::canDownload() const
{
    if ( M_url )
        return true;
    else if ( M_github )
        return true;
    else if ( M_girder && M_girder->canDownload() )
        return true;
    return false;
}

bool
RemoteData::canUpload() const
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

void
RemoteData::upload( std::string const& dataPath, std::string const& parentId ) const
{
    if ( M_girder && M_girder->canUpload() )
        M_girder->upload( dataPath, parentId );
}

std::vector<std::pair<std::string,std::string>>
RemoteData::createFolder( std::string const& folderPath, std::string const& parentId ) const
{
    if ( M_girder && M_girder->canUpload() )
        return M_girder->createFolder( folderPath, parentId );
    return {};
}

RemoteData::URL::URL( std::string const& url, WorldComm const& worldComm )
    :
    M_worldComm( worldComm )
{
    std::regex ex("(http|https)://([^/ :]+):?([^/ ]*)(/?[^ #?]*)\\x3f?([^ #]*)#?([^ ]*)");
    std::cmatch what;
    if( regex_match(url.c_str(), what, ex) )
    {
        M_url = url;
        M_protocol = std::string(what[1].first, what[1].second);
        M_domain   = std::string(what[2].first, what[2].second);
        M_port     = std::string(what[3].first, what[3].second);
        M_path     = std::string(what[4].first, what[4].second);
        M_query    = std::string(what[5].first, what[5].second);
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

bool
RemoteData::URL::isValid() const
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
        if ( p.has_filename() && !Feel::filename_is_dot( p )  && !Feel::filename_is_dot_dot( p ) )
            filename = p.filename().string();
        else
            filename = "download";
    }

    fs::path dir = fs::path(_dir);
    bool addSubDirHashURL = false;
    if ( addSubDirHashURL )
    {
        size_t hashURL = std::hash<std::string>()( M_url );
        std::string subdir = (boost::format("%1%")% hashURL).str();
        dir/=fs::path(subdir);
    }
    std::string thefilename = (dir/filename).string();

    if ( !M_worldComm.isMasterRank() )
    {
        M_worldComm.barrier();
        return thefilename;
    }
    std::string url = M_url;
#if defined(FEELPP_HAS_LIBCURL)
    CURL *curl_handle;

    if ( !fs::exists( dir ) )
        fs::create_directories( dir );

    curl_global_init(CURL_GLOBAL_ALL);

    /* init the curl session */
    curl_handle = curl_easy_init();

    /* set URL to get here */
    curl_easy_setopt(curl_handle, CURLOPT_URL, url.c_str());
#if 0
    /* Switch on full protocol/debug output while testing */
    curl_easy_setopt(curl_handle, CURLOPT_VERBOSE, 1L);

    /* disable progress meter, set to 0L to enable and disable debug output */
    curl_easy_setopt(curl_handle, CURLOPT_NOPROGRESS, 1L);
#endif

    /* send all data to this function  */
    curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, write_data);

    /* open the file */
    std::ofstream pagefile/*ofs*/ ( thefilename, std::ios::out|std::ios::binary);
    if(pagefile) {

        /* write the page body to this file handle */
        curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, &pagefile);

        /* get it! */
        curl_easy_perform(curl_handle);

        /* close the file */
        pagefile.close();
    }

    /* cleanup curl stuff */
    curl_easy_cleanup(curl_handle);

    curl_global_cleanup();

#else
    CHECK( false ) << "LIBCURL is not detected";
#endif
    M_worldComm.barrier();
    return thefilename;
}


pt::ptree
convertDescToPropertyTree( std::string const& desc )
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

    // create the property tree
    pt::ptree pt;
    std::istringstream istr( newDesc );
    pt::read_json( istr, pt );
    return pt;
}


RemoteData::Github::Github( std::string const& desc, WorldComm const& worldComm )
    :
    M_worldComm( worldComm )
{
    std::regex ex("([ ]*)github([ ]*):([ ]*)([{])([^]*)([}])");
    std::cmatch what;
    if( !regex_match(desc.c_str(), what, ex) )
        return;

    CHECK( what.size() == 7 ) << "invalid size";
    std::vector<std::string> keysvalues;
    std::string exprtosplit = std::string(what[5].first, what[5].second);
#if 0
    boost::split( keysvalues, exprtosplit, boost::is_any_of(","), boost::token_compress_on );
    for ( std::string const& keyvalue : keysvalues )
    {
        std::vector<std::string> keyvalueSplitted;
        boost::split( keyvalueSplitted, keyvalue, boost::is_any_of(":"), boost::token_compress_on );

        CHECK( keyvalueSplitted.size() == 2 ) << "invalid size";
        boost::trim(keyvalueSplitted[0]);
        boost::trim(keyvalueSplitted[1]);
        if ( keyvalueSplitted[0] == "owner" ) M_owner = keyvalueSplitted[1];
        else if ( keyvalueSplitted[0] == "repo" ) M_repo = keyvalueSplitted[1];
        else if ( keyvalueSplitted[0] == "branch" ) M_branch = keyvalueSplitted[1];
        else if ( keyvalueSplitted[0] == "path" ) M_path = keyvalueSplitted[1];
        else if ( keyvalueSplitted[0] == "token" ) M_token = keyvalueSplitted[1];
    }
#else
    pt::ptree pt = convertDescToPropertyTree( exprtosplit );
    if ( auto it = pt.get_optional<std::string>("owner") )
        M_owner = *it;
    if ( auto it = pt.get_optional<std::string>("repo") )
        M_repo = *it;
    if ( auto it = pt.get_optional<std::string>("branch") )
        M_branch = *it;
    if ( auto it = pt.get_optional<std::string>("path") )
        M_path = fs::path(*it).remove_trailing_separator().remove_trailing_separator().string();
    if ( auto it = pt.get_optional<std::string>("token") )
        M_token = *it;
    // if token is empty, try looking for it in environment variable FEELPP_GITHUB_TOKEN
    if ( M_token.empty() )
    {
        char* env;
        env = getenv( "FEELPP_GITHUB_TOKEN" );
        if ( env != NULL && env[0] != '\0' )
        {
            M_token = env;
        }
    }
#endif

    if ( M_owner.empty() )
        M_owner = "feelpp";
    if ( M_repo.empty() )
        M_repo = "feelpp";

#if 0
    std::cout << "owner: " << M_owner << "\n"
              << "repo: " << M_repo << "\n"
              << "branch: " << M_branch << "\n"
              << "path: " << M_path << "\n";
#endif
}
bool
RemoteData::Github::isInit() const
{
    return !M_owner.empty() && !M_repo.empty();
}
std::vector<std::string>
RemoteData::Github::download( std::string const& dir ) const
{
    std::vector<std::string> downloadFileOrFolder;
    if ( M_worldComm.isMasterRank() )
    {
        std::string url = "https://api.github.com/repos/" + M_owner + "/" + M_repo +"/contents/" + M_path;
        if ( !M_branch.empty() )
            url += "?ref=" + M_branch;
        std::vector<std::string> headers;
        headers.push_back("Accept: application/vnd.github.v3.json");
        headers.push_back("User-Agent: feelpp-agent");
        if ( !M_token.empty() )
            headers.push_back( "Authorization: token "+M_token );

        std::ostringstream omemfile;
        requestHTTPGET( url,headers, omemfile );
        std::istringstream istr( omemfile.str() );
        pt::ptree ptree;
        pt::read_json(istr, ptree);

        if ( !fs::exists( dir ) )
            fs::create_directories( dir );

        if ( auto testFile = ptree.get_optional<std::string>("path") )
        {
            std::string filename = ptree.get<std::string>("name");
            headers[0] = "Accept: application/vnd.github.v3.raw";
            std::string filepath = (fs::path(dir)/filename).string();
            std::ofstream ofile( filepath, std::ios::out|std::ios::binary);
            requestHTTPGET( url,headers,ofile );
            ofile.close();
            downloadFileOrFolder.push_back( filepath );
        }
        else
        {
            std::string subdir = (M_path.empty())? M_repo : fs::path(M_path).filename().string();
            std::string newdir = (fs::path(dir)/subdir).string();
            fs::create_directories( newdir );

            this->downloadFolderRecursively( ptree, newdir );
            downloadFileOrFolder.push_back( newdir );
        }
    }
    mpi::broadcast( M_worldComm.globalComm(), downloadFileOrFolder, M_worldComm.masterRank() );
    return downloadFileOrFolder;
}

void
RemoteData::Github::downloadFolderRecursively( pt::ptree const& ptree, std::string const& dir ) const
{
    std::vector<std::string> headers;
    headers.push_back("Accept: application/vnd.github.v3.json");
    headers.push_back("User-Agent: feelpp-agent");
    if ( !M_token.empty() )
        headers.push_back( "Authorization: token "+M_token );

    for (auto const& item : ptree)
    {
        std::string type = item.second.get<std::string>("type");
        std::string name = item.second.get<std::string>("name");
        std::string pathInUrl = item.second.get<std::string>("path");
        std::string url = "https://api.github.com/repos/" + M_owner + "/" + M_repo +"/contents/" + pathInUrl;
        if ( !M_branch.empty() )
            url += "?ref=" + M_branch;

        if ( type == "file" )
        {
            // download file
            headers[0] = "Accept: application/vnd.github.v3.raw";
            std::string filepath = (fs::path(dir)/name).string();
            std::ofstream ofile( filepath, std::ios::out|std::ios::binary);
            requestHTTPGET( url,headers,ofile );
            ofile.close();
        }
        else if ( type == "dir" )
        {
            // get ptree subdir
            headers[0] = "Accept: application/vnd.github.v3.json";
            std::ostringstream omemfile;
            requestHTTPGET( url,headers,omemfile );
            std::istringstream istr( omemfile.str() );
            pt::ptree ptreeSubdir;
            pt::read_json( istr,ptreeSubdir );
            // create subdir
            std::string newdir = (fs::path(dir)/name).string();
            fs::create_directories( newdir );
            // recursive call
            this->downloadFolderRecursively( ptreeSubdir, newdir );
        }
    }
}


RemoteData::Girder::Girder( std::string const& desc, WorldComm const& worldComm )
    :
    M_worldComm( worldComm )
{
    std::regex ex("([ ]*)girder([ ]*):([ ]*)([{])([^]*)([}])");
    std::cmatch what;
    if( !regex_match(desc.c_str(), what, ex) )
        return;

    CHECK( what.size() == 7 ) << "invalid size";
    std::vector<std::string> keysvalues;
    std::string exprtosplit = std::string(what[5].first, what[5].second);

    pt::ptree pt = convertDescToPropertyTree( exprtosplit );

    if ( auto it = pt.get_optional<std::string>("url") )
        M_url = *it;
    if ( auto it = pt.get_optional<std::string>("api_key") )
        M_apiKey = *it;
    if ( auto it = pt.get_optional<std::string>("token") )
        M_token = *it;
    if ( auto it = pt.get_child_optional("file") )
    {
        for( auto const& item : pt.get_child("file") )
            M_fileIds.insert( item.second.get_value<std::string>() );
        if ( M_fileIds.empty() )
            M_fileIds.insert( pt.get<std::string>("file") );
    }
    if ( auto it = pt.get_child_optional("folder") )
    {
        for( auto const& item : pt.get_child("folder") )
            M_folderIds.insert( item.second.get_value<std::string>() );
        if ( M_folderIds.empty() )
            M_folderIds.insert( pt.get<std::string>("folder") );
    }

    if ( M_url.empty() )
        M_url = "https://girder.math.unistra.fr";

    if ( M_apiKey.empty() )
    {
        char* env;
        env = getenv( "FEELPP_GIRDER_API_KEY" );
        if ( env != NULL && env[0] != '\0' )
        {
            M_apiKey = env;
        }
    }

    if ( M_token.empty() )
    {
        char* env;
        env = getenv( "FEELPP_GIRDER_TOKEN" );
        if ( env != NULL && env[0] != '\0' )
        {
            M_token = env;
        }
    }

#if 0
    std::cout << "url: " << M_url << "\n";
    for ( std::string const& fileId : M_fileIds )
        std::cout << "file id: " << fileId << "\n";
#endif

}

void
RemoteData::Girder::setFolderIds( std::string const& folderId )
{
    M_folderIds.clear();
    M_folderIds.insert( folderId );
}

bool
RemoteData::Girder::isInit() const
{
    return !M_url.empty();
}
bool
RemoteData::Girder::canDownload() const
{
    return this->isInit() && (!M_fileIds.empty() || !M_folderIds.empty());
}
bool
RemoteData::Girder::canUpload() const
{
    return this->isInit() && ( !M_token.empty() || !M_apiKey.empty() );
}

std::vector<std::string>
RemoteData::Girder::download( std::string const& dir ) const
{
    std::vector<std::string> downloadedFileOrFolder;
    if ( M_worldComm.isMasterRank() )
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
        // delete token if created
        if ( M_token.empty() && !M_apiKey.empty() && !token.empty() )
            this->removeToken( token );
    }
    mpi::broadcast( M_worldComm.globalComm(), downloadedFileOrFolder, M_worldComm.masterRank() );
    return downloadedFileOrFolder;
}

std::string
RemoteData::Girder::downloadFile( std::string const& fileId, std::string const& dir, std::string const& token ) const
{
    std::string downloadedFile;
    // get metadata info
    std::string urlFileInfo = M_url+"/api/v1/file/" + fileId;
    std::vector<std::string> headersFileInfo;
    headersFileInfo.push_back("Accept: application/json");
    if ( !token.empty() )
        headersFileInfo.push_back( "Girder-Token: "+token );
    std::ostringstream omemfile;
    StatusRequestHTTP status = requestHTTPGET( urlFileInfo, headersFileInfo, omemfile );
    if ( !status.success() )
    {
        std::cout << "requestHTTPGET error : " << status.msg() << "\n";
        return {};
    }
    if ( status.code() != 200 )
    {
        if ( status.code() == 400 )
            std::cout << "invalid file id : " << fileId << "\n";
        else if ( status.code() == 403 )
            std::cout << "Read access was denied on the file\n";
        else
            std::cout << "Error with getting metadata (before download)\n";
        return {};
    }
    // convert to property tree
    std::istringstream istr( omemfile.str() );
    pt::ptree pt;
    pt::read_json(istr, pt);

    // extract info of ptree
    auto itFileName = pt.get_optional<std::string>("name");
    CHECK( itFileName ) << "invalid id : not a file or not exists";
    std::string filename = *itFileName;
    auto itMimeType = pt.get_optional<std::string>("mimeType");
    std::string sha512;
    if ( auto itSha512 = pt.get_optional<std::string>("sha512") )
        sha512 = *itSha512;

    std::string filepath = (fs::path(dir)/filename).string();
    std::string metadatapath = (fs::path(dir)/(filename+".metadata.json")).string();

    // download can be ignored if file exists and have same sha512
    bool doDownload = true;
    if ( !sha512.empty() && fs::exists( filepath ) && fs::is_regular_file( filepath ) && fs::exists( metadatapath ) )
    {
        pt::ptree ptmdExisting;
        pt::read_json( metadatapath, ptmdExisting );
        if ( auto itSha512 = ptmdExisting.get_optional<std::string>("sha512") )
        {
            if ( sha512 == *itSha512 )
                doDownload = false;
        }
    }
    // download the file
    if ( doDownload )
    {
        std::string urlFileDownload = M_url+"/api/v1/file/" + fileId + "/download";
        std::vector<std::string> headersFileDownload;
        if ( itMimeType )
            headersFileDownload.push_back( "Accept: " + *itMimeType );
        if ( !token.empty() )
            headersFileDownload.push_back( "Girder-Token: "+token );
        std::ofstream ofile( filepath, std::ios::out|std::ios::binary);
        status = requestHTTPGET( urlFileDownload,headersFileDownload,ofile );
        ofile.close();
        if ( !status.success() )
        {
            std::cout << "requestHTTPGET error : " << status.msg() << "\n";
            return {};
        }
        if ( status.code() != 200 )
        {
            if ( status.code() == 400 )
                std::cout << "invalid file id : " << fileId << "\n";
            else if ( status.code() == 403 )
                std::cout << "Read access was denied on the file\n";
            else
                std::cout << "Error with downloading\n";
            return {};
        }
        // save metadata
        std::ofstream ofileMetadata( metadatapath, std::ios::out);
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
    // get metadata info
    std::string urlFolderInfo = M_url+"/api/v1/folder/" + folderId;
    std::vector<std::string> headersFolderInfo;
    headersFolderInfo.push_back("Accept: application/json");
    if ( !token.empty() )
        headersFolderInfo.push_back( "Girder-Token: "+token );
    std::ostringstream omemfile;
    requestHTTPGET( urlFolderInfo, headersFolderInfo, omemfile );
    // convert to property tree
    std::istringstream istr( omemfile.str() );
    pt::ptree pt;
    pt::read_json(istr, pt);

    // extract info of ptree
    auto itFolderName = pt.get_optional<std::string>("name");
    CHECK( itFolderName ) << "invalid id : not a folder or not exists";
    std::string foldername = *itFolderName;

    // download the folder
    std::string urlFolderDownload = M_url+"/api/v1/folder/" + folderId + "/download";
    std::vector<std::string> headersFolderDownload;
    if ( !token.empty() )
        headersFolderDownload.push_back( "Girder-Token: "+token );
    std::string filepath = (fs::path(dir)/(foldername+".zip")).string();
    std::ofstream ofile( filepath, std::ios::out|std::ios::binary);
    requestHTTPGET( urlFolderDownload,headersFolderDownload,ofile );
    ofile.close();
    // save metadata
    std::string metadatapath = (fs::path(dir)/(foldername+".metadata.json")).string();
    std::ofstream ofileMetadata( metadatapath, std::ios::out);
    ofileMetadata << omemfile.str();
    ofileMetadata.close();

    downloadedFolder = filepath;
    return downloadedFolder;
}

void
RemoteData::Girder::upload( std::string const& dataPath, std::string const& parentId ) const
{
    CHECK( !M_token.empty() || !M_apiKey.empty() ) << "authentication unavailable";
    std::string currentParentId = parentId;
    if ( parentId.empty() && !M_folderIds.empty() )
        currentParentId = *M_folderIds.begin();
    CHECK( !currentParentId.empty() ) << "a parentId is required";

    if ( M_worldComm.isMasterRank() )
    {
        std::string token = M_token;
        if ( M_token.empty() )
            token = this->createToken();

        fs::path dataFsPath( dataPath );
        if( fs::is_directory( dataFsPath ) && Feel::filename_is_dot( dataFsPath.filename() ) )
        {
            fs::directory_iterator end_itr;
            for ( fs::directory_iterator itr( dataFsPath ); itr != end_itr; ++itr )
                this->uploadRecursively( itr->path().string(), currentParentId, token );
        }
        else
            this->uploadRecursively( dataPath, currentParentId, token );

        if ( M_token.empty() && !token.empty() )
            this->removeToken( token );
    }
    M_worldComm.barrier();
}

std::vector<std::pair<std::string,std::string> >
RemoteData::Girder::createFolder( std::string const& folderPath, std::string const& parentId ) const
{
    CHECK( !M_token.empty() || !M_apiKey.empty() ) << "authentication unavailable";
    std::string currentParentId = parentId;
    if ( parentId.empty() && !M_folderIds.empty() )
        currentParentId = *M_folderIds.begin();
    CHECK( !currentParentId.empty() ) << "a parentId is required";

    std::vector<boost::tuple<std::string,std::string> > foldersInfo;
    if ( M_worldComm.isMasterRank() )
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
    mpi::broadcast( M_worldComm.globalComm(), foldersInfo, M_worldComm.masterRank() );
    std::vector<std::pair<std::string,std::string> > res;
    for ( auto const& folderInfo : foldersInfo )
        res.push_back( std::make_pair( boost::get<0>( folderInfo ), boost::get<1>( folderInfo ) ) );
    return res;
}


void
RemoteData::Girder::uploadRecursively( std::string const& dataPath, std::string const& parentId, std::string const& token ) const
{
    fs::path dataFsPath( dataPath );
    if( fs::is_directory( dataFsPath ) )
    {
        std::string folderName = dataFsPath.filename().string();
        std::string newParentId = this->createFolderImpl( folderName, parentId, token );

        fs::directory_iterator end_itr;
        for ( fs::directory_iterator itr( dataFsPath ); itr != end_itr; ++itr )
            this->uploadRecursively( itr->path().string(), newParentId, token );
    }
    else if ( fs::is_regular_file( dataFsPath ) )
    {
        this->uploadFile( dataPath, parentId, token );
    }
}

void
RemoteData::Girder::uploadFile( std::string const& filepath, std::string const& parentId, std::string const& token ) const
{
    std::string filename = fs::path(filepath).filename().string();
    std::string fileExtension = fs::path(filepath).extension().string();

    std::ifstream imemfile( filepath, std::ios::binary );
    // get length of file
    imemfile.seekg( 0, imemfile.end );
    long fsize = imemfile.tellg();
    imemfile.seekg( 0, imemfile.beg );

    std::string urlFileUpload = M_url+"/api/v1/file?parentType=folder&parentId=" + parentId;
    urlFileUpload += "&name="+filename;
    urlFileUpload += "&size="+(boost::format("%1%")%fsize).str();
    std::string mimeType;
    auto itFindMimeType = mimeTypes.find( fileExtension );
    if ( itFindMimeType != mimeTypes.end() )
    {
        std::vector<std::string> exprSplitted;
        boost::split( exprSplitted, itFindMimeType->second, boost::is_any_of("/"), boost::token_compress_on );
        mimeType = exprSplitted[0];
        for ( int k=1;k<exprSplitted.size();++k )
            mimeType += "%2F" + exprSplitted[k];
    }
    if ( !mimeType.empty() )
        urlFileUpload += "&mimeType="+mimeType;

    std::vector<std::string> headersFileInfo;
    headersFileInfo.push_back("Accept: application/json");
    headersFileInfo.push_back("Content-Type: application/form-data");
    CHECK( !token.empty() ) << "a token is required for upload";
    headersFileInfo.push_back( "Girder-Token: "+token );

    std::ostringstream omemfile;
    requestHTTPPOST( urlFileUpload, headersFileInfo, imemfile, fsize, omemfile );
}

std::string
RemoteData::Girder::createFolderImpl( std::string const& folderName, std::string const& parentId, std::string const& token ) const
{
    std::string urlCreateFolder = M_url+"/api/v1/folder?parentType=folder&parentId=" + parentId;
    urlCreateFolder += "&name="+ folderName;
    urlCreateFolder += "&reuseExisting=true";

    std::vector<std::string> headers;
    headers.push_back("Accept: application/json");
    headers.push_back("Content-Type: application/x-www-form-urlencoded");
    CHECK( !token.empty() ) << "a token is required for upload";
    headers.push_back( "Girder-Token: "+token );

    std::ostringstream omemfile;
    requestHTTPPOST( urlCreateFolder, headers, omemfile );
    // convert to property tree
    std::istringstream istr( omemfile.str() );
    pt::ptree pt;
    pt::read_json(istr, pt);

    std::string folderIdCreated;
    if ( auto it = pt.get_optional<std::string>("_id") )
        folderIdCreated = *it;
    return folderIdCreated;
}

std::string
RemoteData::Girder::createToken( int duration ) const
{
    std::string urlCreateToken = M_url+"/api/v1/api_key/token";
    urlCreateToken += "?key=" + M_apiKey;
    urlCreateToken += "&duration=" + std::to_string( duration );

    std::vector<std::string> headers;
    headers.push_back("Accept: application/json");
    headers.push_back("Content-Type: application/json");

    std::ostringstream omemfile;
    StatusRequestHTTP status = requestHTTPPOST( urlCreateToken, headers, omemfile );
    if ( !status.success() )
    {
        std::cout << "requestHTTPGET error : " << status.msg() << "\n";
        return {};
    }
    if ( status.code() != 200 )
    {
        if ( status.code() == 400 )
            std::cout << "invalid api key\n";
        else
            std::cout << "Error with token creation\n";
        return {};
    }

    // convert to property tree
    std::istringstream istr( omemfile.str() );
    pt::ptree pt;
    pt::read_json(istr, pt);

    std::string tokenCreated;
    if ( auto it = pt.get_child_optional( "authToken" ) )
    {
        if ( auto it2 = it->get_optional<std::string>("token") )
            tokenCreated = *it2;
    }
    return tokenCreated;
}

void
RemoteData::Girder::removeToken( std::string const& token ) const
{
    std::string urlRemoveToken = M_url+"/api/v1/token/session";
    std::vector<std::string> headers;
    headers.push_back("Accept: application/json");
    CHECK( !token.empty() ) << "a token is required for upload";
    headers.push_back( "Girder-Token: "+token );

    std::ostringstream omemfile;
    requestHTTPDELETE( urlRemoveToken, headers, omemfile );
}


} // namespace Feel
