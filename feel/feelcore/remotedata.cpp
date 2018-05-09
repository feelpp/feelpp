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

extern "C" {

#include <curl/curl.h>

}
namespace Feel
{

RemoteData::RemoteData( std::string const& desc, WorldComm const& worldComm )
{
    RemoteData::URL urlTool( desc,worldComm );
    if ( urlTool.isValid() )
        M_url.emplace( urlTool );
        //M_url = boost::optional<URL>( urlTool );
}

bool
RemoteData::canDownload() const
{
    if ( M_url )
        return true;
    return false;
}
std::string
RemoteData::download( std::string const& dir, std::string const& filename ) const
{
    if ( M_url )
        return M_url->download( dir, filename );
    return std::string("");
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
        std::cout << "[" << url << "]" << std::endl;
        std::cout << "protocol: " << M_protocol << std::endl;
        std::cout << "domain: " << M_domain << std::endl;
        std::cout << "port:" << M_port << std::endl;
        std::cout << "path:" << M_path << std::endl;
        std::cout << "query:" << M_query << std::endl;
        //std::cout << "isDIr:" << fs::is_directory(M_path) << std::endl;
        std::cout << "filename: " << fs::path(M_path).filename().string() << std::endl;

        size_t hashURL = std::hash<std::string>()( M_url );
        std::cout << "hash:" << hashURL << std::endl;
    }
}

bool
RemoteData::URL::isValid() const
{
    return !M_url.empty();
}





static size_t write_data(char/*void*/ *ptr, size_t size, size_t nmemb, void *stream)
{
    //size_t written = fwrite(ptr, size, nmemb, (FILE *)stream);
    ((std::ofstream * ) stream)->write(ptr,nmemb);
    return nmemb/*written*/;
}


std::string
RemoteData::URL::download( std::string const& _dir, std::string const& _filename ) const
{
    CHECK( this->isValid() ) << "url is not valid";

    std::string filename = _filename;
    if ( filename.empty() )
    {
        fs::path p = fs::path( M_path );
        if ( p.has_filename() && !p.filename_is_dot() && !p.filename_is_dot_dot() )
            filename = p.filename().string();
        else
            filename = "download";
    }
    size_t hashURL = std::hash<std::string>()( M_url );
    std::string subdir = (boost::format("%1%")% hashURL).str();

    fs::path dir = fs::path(_dir)/fs::path(subdir);
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

void
RemoteData::Girder::download( size_type id, std::string const& dir )
{
#if defined(FEELPP_HAS_LIBCURL)
    CURL *curl_handle;
#else
    CHECK( false ) << "LIBCURL is not detected";
#endif
}

} // namespace Feel
