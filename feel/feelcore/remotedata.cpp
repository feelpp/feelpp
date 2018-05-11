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

static size_t write_data(char/*void*/ *ptr, size_t size, size_t nmemb, void *stream)
{
    //size_t written = fwrite(ptr, size, nmemb, (FILE *)stream);
    //((std::ofstream * ) stream)->write(ptr,nmemb);
    ((std::ostream * ) stream)->write(ptr,nmemb);
    return nmemb/*written*/;
}

void requestHTTPGET( std::string const& url, std::vector<std::string> const& headers, std::string& response)
{
#if defined(FEELPP_HAS_LIBCURL)
    curl_global_init(CURL_GLOBAL_ALL);
    CURL *curl_handle;
    curl_handle = curl_easy_init();
    curl_easy_setopt(curl_handle, CURLOPT_URL, url.c_str() );
    curl_easy_setopt(curl_handle, CURLOPT_HTTPGET, 1L);

    struct curl_slist *list = NULL;
    //list = curl_slist_append(list, "Accept: application/vnd.github.v3.raw");
    for ( std::string const& header : headers )
        list = curl_slist_append(list, header.c_str() );
    curl_easy_setopt(curl_handle, CURLOPT_HTTPHEADER, list);

    /* send all data to this function  */
    curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, write_data);

    std::ostringstream ofile;

    curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, &ofile);

    /* get it! */
    curl_easy_perform(curl_handle);

    /* cleanup curl stuff */
    curl_easy_cleanup(curl_handle);

    curl_global_cleanup();

    response = ofile.str();
#else
    CHECK( false ) << "LIBCURL is not detected";
#endif
}



RemoteData::RemoteData( std::string const& desc, WorldComm const& worldComm )
{
    RemoteData::URL urlTool( desc,worldComm );
    if ( urlTool.isValid() )
        M_url.emplace( urlTool );
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

RemoteData::Github::Github( std::string const& desc, WorldComm const& worldComm )
    :
    M_worldComm( worldComm )
{
    std::string subex = "(owner|repo|branch|path):([^]*)";
    std::string subsep = "([ ]*)([,])([ ]*)";
    std::regex ex("([ ]*)github([ ]*):([ ]*)[{]([ ]*)"+subex+subsep+subex+subsep+subex+subsep+subex+"([ ]*)([}])");
    std::cmatch what;
    if( !regex_match(desc.c_str(), what, ex) )
        return;

    CHECK( what.size() == 24 ) << "invalid size";
    for ( int k : { 5,10,15,20 } )
    {
        if ( what[k] == "owner" ) M_owner = what[k+1];
        else if ( what[k] == "repo" ) M_repo = what[k+1];
        else if ( what[k] == "branch" ) M_branch = what[k+1];
        else if ( what[k] == "path" ) M_path = what[k+1];
    }

    std::cout << "owner: " << M_owner << "\n"
              << "repo: " << M_repo << "\n"
              << "branch: " << M_branch << "\n"
              << "path: " << M_path << "\n";
}
void
RemoteData::Github::download( std::string const& dir )
{
    std::string url = "https://api.github.com/repos/" + M_owner + "/" + M_repo +"/contents/" + M_path + "?ref="+M_branch;
    std::vector<std::string> headers;
    headers.push_back("Accept: application/vnd.github.v3.json");
    headers.push_back("User-Agent: feelpp-agent");

    std::string response;
    requestHTTPGET( url,headers, response );
    std::istringstream istr( response );
    pt::ptree ptree;
    pt::read_json(istr, ptree);

    if ( auto testFile = ptree.get_optional<std::string>("path") )
    {
        // a file
        std::string filename = ptree.get<std::string>("name");
        std::string fileurl = ptree.get<std::string>("download_url");
        RemoteData::URL urlTool(fileurl,M_worldComm);
        urlTool.download( dir,filename );
    }
    else
    {
        // a folder
        for (auto& item : ptree)
            std::cout << "THE PATH : " << item.second.get<std::string>("path") << "\n";
    }
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
