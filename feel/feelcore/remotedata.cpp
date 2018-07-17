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

StatusRequestHTTP requestHTTPPOST( std::string const& url, std::vector<std::string> const& headers, std::istream & ifile, long fsize, std::ostream & ofile )
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

    // post file
#if 1
    char * postthis = new char [fsize];
    ifile.read( postthis,fsize );
    res = curl_easy_setopt(curl_handle, CURLOPT_POSTFIELDS, postthis);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );
#else
    // very slow (need to understand)
    res = curl_easy_setopt(curl_handle, CURLOPT_POST, 1L);
    res = curl_easy_setopt(curl_handle, CURLOPT_READDATA, &ifile);
    res = curl_easy_setopt(curl_handle, CURLOPT_READFUNCTION, read_data);
#endif
    res = curl_easy_setopt(curl_handle, CURLOPT_POSTFIELDSIZE, fsize);
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

    delete [] postthis;
    return StatusRequestHTTP( true, http_code, "" );
#else
    CHECK( false ) << "LIBCURL is not detected";
    return StatusRequestHTTP( false );
#endif
}

StatusRequestHTTP requestHTTPCUSTOM( std::string const& customRequest, std::string const& url, std::vector<std::string> const& headers, std::ostream & ofile )
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
    // request type
    res = curl_easy_setopt(curl_handle, CURLOPT_CUSTOMREQUEST, customRequest.c_str()/*"DELETE"*/);
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

StatusRequestHTTP requestDownloadURL( std::string const& url, std::ostream & ofile)
{
#if defined(FEELPP_HAS_LIBCURL)
    CURLcode res;
    res = curl_global_init(CURL_GLOBAL_ALL);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );

    /* init the curl session */
    CURL *curl_handle = curl_easy_init();
    if ( !curl_handle ) return StatusRequestHTTP( false, "fail to run curl_easy_init" );

    /* set URL to get here */
    res = curl_easy_setopt(curl_handle, CURLOPT_URL, url.c_str());
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );
#if 0
    /* Switch on full protocol/debug output while testing */
    res = curl_easy_setopt(curl_handle, CURLOPT_VERBOSE, 1L);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );
    /* disable progress meter, set to 0L to enable and disable debug output */
    res = curl_easy_setopt(curl_handle, CURLOPT_NOPROGRESS, 1L);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );
#endif

    /* send all data to this function  */
    res = curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, write_data);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );

    /* write the page body to this file handle */
    res = curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, &ofile);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );

    /* get it! */
    res = curl_easy_perform(curl_handle);
    if ( res != CURLE_OK ) return StatusRequestHTTP( false, curl_easy_strerror( res ) );

    long http_code = 0;
    res = curl_easy_getinfo(curl_handle, CURLINFO_RESPONSE_CODE, &http_code);

    /* cleanup curl stuff */
    curl_easy_cleanup(curl_handle);
    curl_global_cleanup();

    return StatusRequestHTTP( true, http_code, "" );
#else
    CHECK( false ) << "LIBCURL is not detected";
    return StatusRequestHTTP( false );
#endif

}


RemoteData::RemoteData( std::string const& desc, WorldComm const& worldComm )
    :
    M_worldComm( worldComm )
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

std::vector<std::string>
RemoteData::upload( std::string const& dataPath, std::string const& parentId, bool sync ) const
{
    if ( M_girder && M_girder->canUpload() )
        return M_girder->upload( dataPath, parentId, sync );
    return {};
}

std::vector<std::vector<std::string>>
RemoteData::upload( std::vector<std::pair<std::string,std::string> > const& dataToUpload, bool sync ) const
{
    if ( M_girder && M_girder->canUpload() )
        return M_girder->upload( dataToUpload, sync );
    return {};
}

void
RemoteData::replaceFile( std::string const& filePath, std::string const& fileId ) const
{
    if ( M_girder && M_girder->canUpload() )
        return M_girder->replaceFile( filePath, fileId );
}

void
RemoteData::replaceFile( std::vector<std::pair<std::string,std::string> > const& filesToReplace ) const
{
    if ( M_girder && M_girder->canUpload() )
        return M_girder->replaceFile( filesToReplace );
}

std::vector<std::pair<std::string,std::string>>
RemoteData::createFolder( std::string const& folderPath, std::string const& parentId, bool sync ) const
{
    if ( M_girder && M_girder->canUpload() )
        return M_girder->createFolder( folderPath, parentId, sync );
    return {};
}

std::tuple<std::vector<std::shared_ptr<RemoteData::FolderInfo>>,std::vector<std::shared_ptr<RemoteData::ItemInfo>>,std::vector<std::shared_ptr<RemoteData::FileInfo>>>
RemoteData::contents() const
{
    if ( M_girder && M_girder->isInit() )
        return M_girder->contents();
    return std::make_tuple(std::vector<std::shared_ptr<FolderInfo>>(),std::vector<std::shared_ptr<ItemInfo>>(),std::vector<std::shared_ptr<FileInfo>>());
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

    std::string url = M_url;

    if ( M_worldComm.isMasterRank() )
    {
        if ( !fs::exists( dir ) )
            fs::create_directories( dir );

        std::ofstream ofile( thefilename, std::ios::out|std::ios::binary);
        /* open the file */
        if( ofile )
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

    M_worldComm.barrier();
    return thefilename;
}


std::pair<bool,pt::ptree>
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
    try {
        pt::read_json( istr, pt );
    }
    catch ( pt::ptree_error const& e )
    {
        return std::make_pair( false, pt );
    }
    return std::make_pair( true, pt );
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
    auto resConvertion = convertDescToPropertyTree( exprtosplit );
    if ( !resConvertion.first )
    {
        if ( M_worldComm.isMasterRank() )
            std::cout << "Github desc has a syntax error : " << desc << "\n";
        return;
    }
    pt::ptree pt = resConvertion.second;
    if ( auto it = pt.get_optional<std::string>("owner") )
        M_owner = *it;
    if ( auto it = pt.get_optional<std::string>("repo") )
        M_repo = *it;
    if ( auto it = pt.get_optional<std::string>("branch") )
        M_branch = *it;
    if ( auto it = pt.get_optional<std::string>("path") )
    {
        if ( Feel::filename_is_dot( fs::path( *it ).filename() ) )
            M_path = fs::path( *it ).parent_path().string();
        else
            M_path = *it;
    }
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
        downloadFileOrFolder = this->downloadImpl( dir );
    }
    mpi::broadcast( M_worldComm.globalComm(), downloadFileOrFolder, M_worldComm.masterRank() );
    return downloadFileOrFolder;
}

std::vector<std::string>
RemoteData::Github::downloadImpl( std::string const& dir ) const
{
    std::vector<std::string> downloadFileOrFolder;

    std::string url = "https://api.github.com/repos/" + M_owner + "/" + M_repo +"/contents/" + M_path;
    if ( !M_branch.empty() )
        url += "?ref=" + M_branch;
    std::vector<std::string> headers;
    headers.push_back("Accept: application/vnd.github.v3.json");
    headers.push_back("User-Agent: feelpp-agent");
    if ( !M_token.empty() )
        headers.push_back( "Authorization: token "+M_token );

    std::ostringstream omemfile;
    StatusRequestHTTP status = requestHTTPGET( url,headers, omemfile );
    if ( !status.success() )
    {
        std::cout << "Github error in requestHTTPGET : " << status.msg() << "\n";
        return {};
    }

    std::istringstream istr( omemfile.str() );
    pt::ptree ptree;
    pt::read_json(istr, ptree);

    if ( status.code() != 200 )
    {
        std::cout << Github::errorMessage( ptree,"get metadata fails", status.code() ) << "\n";
        return {};
    }

    if ( !fs::exists( dir ) )
        fs::create_directories( dir );

    if ( auto testFile = ptree.get_optional<std::string>("path") )
    {
        std::string filename = ptree.get<std::string>("name");
        headers[0] = "Accept: application/vnd.github.v3.raw";
        std::string filepath = (fs::path(dir)/filename).string();
        std::ofstream ofile( filepath, std::ios::out|std::ios::binary);
        status = requestHTTPGET( url,headers,ofile );
        ofile.close();
        if ( !status.success() )
        {
            std::cout << "Github error in requestHTTPGET : " << status.msg() << "\n";
            return {};
        }
        if ( status.code() != 200 )
        {
            std::cout << Github::errorMessage( pt::ptree{}, "download file fails", status.code() ) << "\n";
            return {};
        }

        downloadFileOrFolder.push_back( filepath );
    }
    else
    {
        std::string subdir = (M_path.empty())? M_repo : fs::path(M_path).filename().string();
        std::string newdir = (fs::path(dir)/subdir).string();
        fs::create_directories( newdir );

        auto resFolder = this->downloadFolderRecursively( ptree, newdir );
        if ( !std::get<0>( resFolder ) )
        {
            std::cout << std::get<1>( resFolder ) << "\n";
            return {};
        }
        downloadFileOrFolder.push_back( newdir );
    }

    return downloadFileOrFolder;
}

std::tuple<bool,std::string>
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
            StatusRequestHTTP status = requestHTTPGET( url,headers,ofile );
            ofile.close();
            if ( !status.success() )
                return std::make_tuple( false, "Github error in requestHTTPGET : " + status.msg() );
            if ( status.code() != 200 )
                return std::make_tuple( false, Github::errorMessage( pt::ptree{},"download file fails", status.code() ) );
        }
        else if ( type == "dir" )
        {
            // get ptree subdir
            headers[0] = "Accept: application/vnd.github.v3.json";
            std::ostringstream omemfile;
            StatusRequestHTTP status = requestHTTPGET( url,headers,omemfile );
            if ( !status.success() )
                return std::make_tuple( false, "Github error in requestHTTPGET : " + status.msg() );
            std::istringstream istr( omemfile.str() );
            pt::ptree ptreeSubdir;
            pt::read_json( istr,ptreeSubdir );
            if ( status.code() != 200 )
                return std::make_tuple( false, Github::errorMessage( ptreeSubdir, "get metadata fails", status.code() ) );
            // create subdir
            std::string newdir = (fs::path(dir)/name).string();
            fs::create_directories( newdir );
            // recursive call
            auto resRecur = this->downloadFolderRecursively( ptreeSubdir, newdir );
            if ( !std::get<0>( resRecur ) )
                return std::make_tuple( false, std::get<1>( resRecur ) );
        }
    }
    return std::make_tuple( true, "" );
}

std::string
RemoteData::Github::errorMessage( pt::ptree const& ptree, std::string const& defaultMsg, uint16_type statusCode )
{
    std::string errMsg = defaultMsg;
    if ( auto itMsg = ptree.get_optional<std::string>("message") )
        errMsg = *itMsg;
    std::string docUrl;
    if ( auto itDoc = ptree.get_optional<std::string>("documentation_url") )
        docUrl = *itDoc;
    std::string errMsgInfo;
    if ( statusCode != invalid_uint16_type_value )
        errMsgInfo += (boost::format("code %1%")%statusCode).str();
    std::string res = "Girder error";
    if ( !errMsgInfo.empty() )
        res += " ("+errMsgInfo+")";
    if ( !errMsg.empty() )
        res += " : " + errMsg;
    if ( !docUrl.empty() )
        res += (boost::format(" [ doc_url: %1% ]")%docUrl).str();
    return res;
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

    auto resConvertion = convertDescToPropertyTree( exprtosplit );
    if ( !resConvertion.first )
    {
        if ( M_worldComm.isMasterRank() )
            std::cout << "Girder desc has a syntax error : " << desc << "\n";
        return;
    }
    pt::ptree pt = resConvertion.second;

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
    if ( auto it = pt.get_child_optional("item") )
    {
        for( auto const& item : pt.get_child("item") )
            M_itemIds.insert( item.second.get_value<std::string>() );
        if ( M_itemIds.empty() )
            M_itemIds.insert( pt.get<std::string>("item") );
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
RemoteData::Girder::errorMessage( pt::ptree const& ptree, std::string const& defaultMsg, uint16_type statusCode )
{
    std::string errMsg = defaultMsg;
    if ( auto itMsg = ptree.get_optional<std::string>("message") )
        errMsg = *itMsg;
    std::string errType;
    if ( auto itType = ptree.get_optional<std::string>("type") )
        errType = *itType;
    std::string errMsgInfo;
    if ( statusCode != invalid_uint16_type_value )
        errMsgInfo += (boost::format("code %1%")%statusCode).str();
    if ( !errMsgInfo.empty() )
        errMsgInfo += ",";
    if ( !errType.empty() )
        errMsgInfo += (boost::format("type %1%")%errType).str();
    std::string res = "Girder error";
    if ( !errMsgInfo.empty() )
        res += " ("+errMsgInfo+")";
    if ( !errMsg.empty() )
        res += " : " + errMsg;
    return res;
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
        std::cout << "Girder error in requestHTTPGET : " << status.msg() << "\n";
        return {};
    }
    // convert to property tree
    std::istringstream istr( omemfile.str() );
    pt::ptree pt;
    pt::read_json(istr, pt);

    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( pt,"getting metadata (before download) fails", status.code() ) << "\n";
        return {};
    }

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
            std::cout << "Girder error in requestHTTPGET : " << status.msg() << "\n";
            return {};
        }
        if ( status.code() != 200 )
        {
            std::cout << Girder::errorMessage( pt::ptree{}, "downloading file fails", status.code() ) << "\n";
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
    StatusRequestHTTP status = requestHTTPGET( urlFolderInfo, headersFolderInfo, omemfile );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPGET : " << status.msg() << "\n";
        return {};
    }

    // convert to property tree
    std::istringstream istr( omemfile.str() );
    pt::ptree pt;
    pt::read_json(istr, pt);

    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( pt,"getting metadata (before download) fails", status.code() ) << "\n";
        return {};
    }

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
    status = requestHTTPGET( urlFolderDownload,headersFolderDownload,ofile );
    ofile.close();
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPGET : " << status.msg() << "\n";
        return {};
    }
    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( pt,"downloading file fails", status.code() ) << "\n";
        return {};
    }
    // save metadata
    std::string metadatapath = (fs::path(dir)/(foldername+".metadata.json")).string();
    std::ofstream ofileMetadata( metadatapath, std::ios::out);
    ofileMetadata << omemfile.str();
    ofileMetadata.close();

    downloadedFolder = filepath;
    return downloadedFolder;
}

std::vector<std::string>
RemoteData::Girder::upload( std::string const& dataPath, std::string const& parentId, bool sync ) const
{
    auto res = this->upload( std::vector<std::pair<std::string,std::string>>(1, std::make_pair(dataPath,parentId) ), sync );
    CHECK( res.size() == 1 ) << "wrong size "<< res.size() << " : must be 1";
    return res[0];
}

std::vector<std::vector<std::string>>
RemoteData::Girder::upload( std::vector<std::pair<std::string,std::string> > const& dataToUpload, bool sync ) const
{
    if ( dataToUpload.empty() )
        return {};

    CHECK( !M_token.empty() || !M_apiKey.empty() ) << "authentication unavailable";

    std::vector<std::vector<std::string>> res( dataToUpload.size() );

    if ( M_worldComm.isMasterRank() )
    {
        // use token if given else create token if api key given
        std::string token = M_token;
        if ( M_token.empty() )
            token = this->createToken();

        for ( int k=0;k<dataToUpload.size();++k )
        {
            std::string dataPath = dataToUpload[k].first;
            std::string parentId = dataToUpload[k].second;
            std::string parentFolderId = parentId;
            std::string parentFileId;
            if ( parentId.empty() )
            {
                if ( !M_folderIds.empty() )
                    parentFolderId = *M_folderIds.begin();
                else if ( !M_fileIds.empty() )
                    parentFileId = *M_fileIds.begin();
            }
            fs::path dataFsPath( dataPath );
            if ( !parentFolderId.empty() ) // upload into a folder
            {
                if( fs::is_directory( dataFsPath ) && Feel::filename_is_dot( dataFsPath.filename() ) )
                {
                    fs::directory_iterator end_itr;
                    for ( fs::directory_iterator itr( dataFsPath ); itr != end_itr; ++itr )
                        res[k] = this->uploadRecursively( itr->path().string(), parentFolderId, token );
                }
                else
                    res[k] = this->uploadRecursively( dataPath, parentFolderId, token );
            }
            else if ( !parentFileId.empty() ) // replace file
            {
                this->replaceFileImpl( dataPath, parentFileId, token );
                res[k].resize(1);
                res[k][0] = parentFileId;
            }
            else
                CHECK( false ) << "a parentId is required";
        }

        // delete token if created
        if ( M_token.empty() && !token.empty() )
            this->removeToken( token );
    }

    if ( sync )
        mpi::broadcast( M_worldComm.globalComm(), res, M_worldComm.masterRank() );

    return res;
}

void
RemoteData::Girder::replaceFile( std::string const& filePath, std::string const& fileId ) const
{
    this->replaceFile( std::vector<std::pair<std::string,std::string>>(1, std::make_pair(filePath,fileId) ) );
}

void
RemoteData::Girder::replaceFile( std::vector<std::pair<std::string,std::string> > const& filesToReplace ) const
{
    if ( filesToReplace.empty() )
        return;

    CHECK( !M_token.empty() || !M_apiKey.empty() ) << "authentication unavailable";

    if ( M_worldComm.isMasterRank() )
    {
        // use token if given else create token if api key given
        std::string token = M_token;
        if ( M_token.empty() )
            token = this->createToken();
        // call replace file request for each file
        for ( int k=0;k<filesToReplace.size();++k )
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

std::vector<std::pair<std::string,std::string> >
RemoteData::Girder::createFolder( std::string const& folderPath, std::string const& parentId, bool sync ) const
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
    if ( sync )
        mpi::broadcast( M_worldComm.globalComm(), foldersInfo, M_worldComm.masterRank() );
    std::vector<std::pair<std::string,std::string> > res;
    for ( auto const& folderInfo : foldersInfo )
        res.push_back( std::make_pair( boost::get<0>( folderInfo ), boost::get<1>( folderInfo ) ) );
    return res;
}


std::vector<std::string>
RemoteData::Girder::uploadRecursively( std::string const& dataPath, std::string const& parentId, std::string const& token ) const
{
    std::vector<std::string> res;
    fs::path dataFsPath( dataPath );
    if( fs::is_directory( dataFsPath ) )
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
        std::string fileUploaded = this->uploadFileImpl( dataPath, parentId, token );
        if ( !fileUploaded.empty() )
            res.push_back( fileUploaded );
    }
    return res;
}

std::string
RemoteData::Girder::uploadFileImpl( std::string const& filepath, std::string const& parentId, std::string const& token ) const
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
    StatusRequestHTTP status = requestHTTPPOST( urlFileUpload, headersFileInfo, imemfile, fsize, omemfile );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPPOST : " << status.msg() << "\n";
        return {};
    }
    // convert to property tree
    std::istringstream istr( omemfile.str() );
    pt::ptree pt;
    pt::read_json(istr, pt);
    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( pt,"upload file fails", status.code() ) << "\n";
        return {};
    }
    std::string fileIdCreated;
    if ( auto it = pt.get_optional<std::string>("_id") )
        fileIdCreated = *it;
    return fileIdCreated;
}

void
RemoteData::Girder::replaceFileImpl( std::string const& filePath, std::string const& fileId, std::string const& token ) const
{
    CHECK( !token.empty() ) << "a token is required for upload";
    CHECK( fs::is_regular_file( filePath ) ) << "must be a file";

    std::ifstream imemfile( filePath, std::ios::binary );
    // get length of file
    imemfile.seekg( 0, imemfile.end );
    long fsize = imemfile.tellg();
    imemfile.seekg( 0, imemfile.beg );

    // request http PUT : contents
    std::string urlFileContents = M_url+"/api/v1/file/"+fileId+"/contents";
    urlFileContents += "?size="+(boost::format("%1%")%fsize).str();
    std::vector<std::string> headersFileContents;
    headersFileContents.push_back("Accept: application/json");
    headersFileContents.push_back("Content-Type: application/json");
    headersFileContents.push_back("Content-Length: 0");
    headersFileContents.push_back( "Girder-Token: "+token );

    std::ostringstream omemfileFileContents;
    StatusRequestHTTP status = requestHTTPCUSTOM( "PUT",urlFileContents, headersFileContents, omemfileFileContents );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPCUSTOM (PUT) : " << status.msg() << "\n";
        return;
    }
    // convert to property tree
    std::istringstream istrFileContents( omemfileFileContents.str() );
    pt::ptree ptFileContents;
    pt::read_json(istrFileContents, ptFileContents);
    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( ptFileContents,"replace file (contents) fails", status.code() ) << "\n";
        return;
    }

    std::string uploadId;
    if ( auto it = ptFileContents.get_optional<std::string>("_id") )
        uploadId = *it;
    CHECK( !uploadId.empty() ) << "no _id in metada returned";


    // request http POST : chunk
    std::string urlFileChunk = M_url+"/api/v1/file/chunk";
    urlFileChunk += "?uploadId="+uploadId;
    urlFileChunk += "&offset=0";
    std::vector<std::string> headersFileChunk;
    headersFileChunk.push_back("Accept: application/json");
    headersFileChunk.push_back("Content-Type: application/form-data");
    headersFileChunk.push_back( "Girder-Token: "+token );

    std::ostringstream omemfile;
    status = requestHTTPPOST( urlFileChunk, headersFileChunk, imemfile, fsize, omemfile );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPPOST : " << status.msg() << "\n";
        return;
    }
    if ( status.code() != 200 )
    {
        // convert to property tree
        std::istringstream istr( omemfile.str() );
        pt::ptree pt;
        pt::read_json(istr, pt);

        std::cout << Girder::errorMessage( pt,"upload file (chunk) fails", status.code() ) << "\n";
        return;
    }

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
    StatusRequestHTTP status = requestHTTPPOST( urlCreateFolder, headers, omemfile );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPPOST : " << status.msg() << "\n";
        return {};
    }

    // convert to property tree
    std::istringstream istr( omemfile.str() );
    pt::ptree pt;
    pt::read_json(istr, pt);
    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( pt,"upload file fails", status.code() ) << "\n";
        return {};
    }

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
        std::cout << "Girder error in requestHTTPPOST : " << status.msg() << "\n";
        return {};
    }

    // convert to property tree
    std::istringstream istr( omemfile.str() );
    pt::ptree pt;
    pt::read_json(istr, pt);

    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( pt,"create token fails", status.code() ) << "\n";
        return {};
    }

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
    StatusRequestHTTP status = requestHTTPCUSTOM( "DELETE",urlRemoveToken, headers, omemfile );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPCUSTOM (DELETE) : " << status.msg() << "\n";
        return;
    }
    if ( status.code() != 200 )
    {
        // convert to property tree
        std::istringstream istr( omemfile.str() );
        pt::ptree pt;
        pt::read_json(istr, pt);

        std::cout << Girder::errorMessage( pt,"remove token fails", status.code() ) << "\n";
        return;
    }
}


std::tuple<std::vector<std::shared_ptr<RemoteData::FolderInfo>>,std::vector<std::shared_ptr<RemoteData::ItemInfo>>,std::vector<std::shared_ptr<RemoteData::FileInfo>>>
RemoteData::Girder::contents() const
{
    auto res = std::make_tuple( std::vector<std::shared_ptr<FolderInfo>>(),
                                std::vector<std::shared_ptr<ItemInfo>>(),
                                std::vector<std::shared_ptr<FileInfo>>() );
    if ( M_worldComm.isMasterRank() )
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
    //if ( sync )

    return res;
}

std::shared_ptr<RemoteData::FileInfo>
RemoteData::Girder::fileInfoImpl( std::string const& fileId, std::string const& token ) const
{
    std::shared_ptr<FileInfo> res;

    std::string url = M_url+"/api/v1/file/" + fileId;

    std::vector<std::string> headers;
    headers.push_back("Accept: application/json");
    if ( !token.empty() )
        headers.push_back( "Girder-Token: "+token );

    std::ostringstream omemfile;
    StatusRequestHTTP status = requestHTTPGET( url, headers, omemfile );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPGET : " << status.msg() << "\n";
        return res;
    }

    // convert to property tree
    std::istringstream istr( omemfile.str() );
    pt::ptree pt;
    pt::read_json(istr, pt);
    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( pt,"file info fails", status.code() ) << "\n";
        return res;
    }

    std::string name, id, mimeType, sha512;
    size_type size = invalid_size_type_value;
    if ( auto it = pt.get_optional<std::string>("_id") )
        id = *it;
    if ( auto it = pt.get_optional<std::string>("name") )
        name = *it;
    if ( auto it = pt.get_optional<size_type>("size") )
        size = *it;
    if ( auto it = pt.get_optional<std::string>("mimeType") )
        mimeType = *it;
    if ( auto it = pt.get_optional<std::string>("sha512") )
        sha512 = *it;
    res = std::make_shared<FileInfo>( name, id, size );
    res->setMimeType( mimeType );
    if ( !sha512.empty() )
        res->setChecksum("sha512",sha512 );
    return res;
}
std::shared_ptr<RemoteData::FolderInfo>
RemoteData::Girder::folderInfoImpl( std::string const& folderId, std::string const& token ) const
{
    std::shared_ptr<FolderInfo> res;

    std::string url = M_url+"/api/v1/folder/" + folderId;

    std::vector<std::string> headers;
    headers.push_back("Accept: application/json");
    if ( !token.empty() )
        headers.push_back( "Girder-Token: "+token );

    std::ostringstream omemfile;
    StatusRequestHTTP status = requestHTTPGET( url, headers, omemfile );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPGET : " << status.msg() << "\n";
        return res;
    }

    // convert to property tree
    std::istringstream istr( omemfile.str() );
    pt::ptree pt;
    pt::read_json(istr, pt);
    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( pt,"folder info fails", status.code() ) << "\n";
        return res;
    }

    std::string name, id;
    size_type size = invalid_size_type_value;
    if ( auto it = pt.get_optional<std::string>("_id") )
        id = *it;
    if ( auto it = pt.get_optional<std::string>("name") )
        name = *it;
    if ( auto it = pt.get_optional<size_type>("size") )
        size = *it;
    res = std::make_shared<FolderInfo>( name, id, size );
    return res;
}

std::shared_ptr<RemoteData::ItemInfo>
RemoteData::Girder::itemInfoImpl( std::string const& itemId, std::string const& token ) const
{
    std::shared_ptr<ItemInfo> res;

    std::string url = M_url+"/api/v1/item/" + itemId;

    std::vector<std::string> headers;
    headers.push_back("Accept: application/json");
    if ( !token.empty() )
        headers.push_back( "Girder-Token: "+token );

    std::ostringstream omemfile;
    StatusRequestHTTP status = requestHTTPGET( url, headers, omemfile );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPGET : " << status.msg() << "\n";
        return res;
    }

    // convert to property tree
    std::istringstream istr( omemfile.str() );
    pt::ptree pt;
    pt::read_json(istr, pt);
    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( pt,"item info fails", status.code() ) << "\n";
        return res;
    }

    std::string name, id;
    size_type size = invalid_size_type_value;
    if ( auto it = pt.get_optional<std::string>("_id") )
        id = *it;
    if ( auto it = pt.get_optional<std::string>("name") )
        name = *it;
    if ( auto it = pt.get_optional<size_type>("size") )
        size = *it;
    res = std::make_shared<ItemInfo>( name, id, size );
    return res;
}

std::shared_ptr<RemoteData::FolderInfo>
RemoteData::Girder::folderContentsImpl( std::string const& folderId, std::string const& token ) const
{
    auto res = folderInfoImpl( folderId, token );
    if ( !res )
        return res;

    std::string urlFolder = M_url+"/api/v1/folder?parentType=folder";
    urlFolder += "&parentId=" + folderId;;
    urlFolder += "&sort=lowerName&sortdir=1";

    std::vector<std::string> headers;
    headers.push_back("Accept: application/json");
    if ( !token.empty() )
        headers.push_back( "Girder-Token: "+token );

    std::ostringstream omemfileFolder;
    StatusRequestHTTP status = requestHTTPGET( urlFolder, headers, omemfileFolder );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPGET : " << status.msg() << "\n";
        return res;
    }

    // convert to property tree
    std::istringstream istrFolder( omemfileFolder.str() );
    pt::ptree ptFolder;
    pt::read_json(istrFolder, ptFolder);
    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( ptFolder,"folder contents (folder) fails", status.code() ) << "\n";
        return res;
    }

    for (auto const& item : ptFolder )
    {
        std::string name, id;
        size_type size = invalid_size_type_value;
        if ( auto it = item.second.get_optional<std::string>("_id") )
            id = *it;
        if ( auto it = item.second.get_optional<std::string>("name") )
            name = *it;
        if ( auto it = item.second.get_optional<size_type>("size") )
            size = *it;
        res->addFolder( std::make_shared<FolderInfo>( name,id,size ) );
    }


    std::string urlItem = M_url+"/api/v1/item";
    urlItem += "?folderId=" + folderId;
    urlItem += "&limit=0&sort=lowerName&sortdir=1";

    std::ostringstream omemfileItem;
    status = requestHTTPGET( urlItem, headers, omemfileItem );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPGET : " << status.msg() << "\n";
        return res;
    }

    // convert to property tree
    std::istringstream istrItem( omemfileItem.str() );
    pt::ptree ptItem;
    pt::read_json(istrItem, ptItem);
    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( ptItem,"folder contents (item) fails", status.code() ) << "\n";
        return res;
    }

    for (auto const& item : ptItem )
    {
        std::string name, id;
        size_type size = invalid_size_type_value;
        if ( auto it = item.second.get_optional<std::string>("_id") )
            id = *it;
        if ( auto it = item.second.get_optional<std::string>("name") )
            name = *it;
        if ( auto it = item.second.get_optional<size_type>("size") )
            size = *it;
        auto itemInfo = std::make_shared<ItemInfo>( name,id,size );
        this->updateFilesImpl( itemInfo, token );
        res->addItem( itemInfo );
    }

    return res;
}

void
RemoteData::Girder::updateFilesImpl( std::shared_ptr<RemoteData::ItemInfo> itemInfo, std::string const& token ) const
{
    if ( !itemInfo )
        return;
    std::string itemId = itemInfo->id();
    if ( itemId.empty() )
        return;

    std::string url = M_url+"/api/v1/item/" + itemId + "/files?limit=0&sort=name&sortdir=1";

    std::vector<std::string> headers;
    headers.push_back("Accept: application/json");
    if ( !token.empty() )
        headers.push_back( "Girder-Token: "+token );

    std::ostringstream omemfile;
    StatusRequestHTTP status = requestHTTPGET( url, headers, omemfile );
    if ( !status.success() )
    {
        std::cout << "Girder error in requestHTTPGET : " << status.msg() << "\n";
        return;
    }

    // convert to property tree
    std::istringstream istr( omemfile.str() );
    pt::ptree pt;
    pt::read_json(istr, pt);
    if ( status.code() != 200 )
    {
        std::cout << Girder::errorMessage( pt,"fileInItem fails", status.code() ) << "\n";
        return;
    }

    for (auto const& ptFile : pt )
    {
        std::string name, id, mimeType, sha512;
        size_type size = invalid_size_type_value;
        if ( auto it = ptFile.second.get_optional<std::string>("_id") )
            id = *it;
        if ( auto it = ptFile.second.get_optional<std::string>("name") )
            name = *it;
        if ( auto it = ptFile.second.get_optional<size_type>("size") )
            size = *it;
        if ( auto it = ptFile.second.get_optional<std::string>("mimeType") )
            mimeType = *it;
        if ( auto it = ptFile.second.get_optional<std::string>("sha512") )
            sha512 = *it;
        auto fileInfo = std::make_shared<FileInfo>( name, id, size );
        fileInfo->setMimeType( mimeType );
        if ( !sha512.empty() )
            fileInfo->setChecksum("sha512",sha512 );
        itemInfo->add( fileInfo );
    }
}

std::ostringstream
RemoteData::FolderInfo::print( size_t nTab ) const
{
    std::ostringstream ostr;
    std::string tab;
    for (int k=0;k<nTab;++k )
        tab += "|   ";
    ostr << tab << "|-- "<< this->name() << "  [type:folder";
    if ( !this->id().empty() )
        ostr << ", id:" << this->id();
    if ( this->size() != invalid_size_type_value )
        ostr << ", size:" << this->size();
    ostr << "]\n";

    for ( auto const& folderInfo : M_folders )
        ostr << folderInfo->print( nTab+1 ).str();
    for ( auto const& itemInfo : M_items )
        ostr << itemInfo->print( nTab+1 ).str();
    for ( auto const& fileInfo : M_files )
        ostr << fileInfo->print( nTab+1 ).str();
    return ostr;
}

std::ostringstream
RemoteData::ItemInfo::print( size_t nTab ) const
{
    std::ostringstream ostr;
    std::string tab;
    for (int k=0;k<nTab;++k )
        tab += "|   ";
    ostr << tab << "|-- "<< this->name() << "  [type:item";
    if ( !this->id().empty() )
        ostr << ", id:" << this->id();
    if ( this->size() != invalid_size_type_value )
        ostr << ", size:" << this->size();
    ostr << "]\n";

    for ( auto const& fileInfo : M_files )
        ostr << fileInfo->print( nTab+1, true ).str();
    return ostr;
}
std::ostringstream
RemoteData::FileInfo::print( size_t nTab, bool isFileInItem ) const
{
    std::ostringstream ostr;
    std::string tab;
    for (int k=0;k<nTab;++k )
        tab += "|   ";
    ostr << tab;
    if ( !isFileInItem )
        ostr << "|-- ";
    else
        ostr << "* ";
    ostr << this->name() << "  [type:file";
    if ( !this->id().empty() )
        ostr << ", id:" << this->id();
    if ( this->size() != invalid_size_type_value )
        ostr << ", size:" << this->size();
    ostr << "]\n";
    return ostr;
}

} // namespace Feel
