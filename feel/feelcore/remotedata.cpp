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

static size_t write_data(char/*void*/ *ptr, size_t size, size_t nmemb, void *stream)
{
    //size_t written = fwrite(ptr, size, nmemb, (FILE *)stream);
    //((std::ofstream * ) stream)->write(ptr,nmemb);
    ((std::ostream * ) stream)->write(ptr,nmemb);
    return nmemb/*written*/;
}

void requestHTTPGET( std::string const& url, std::vector<std::string> const& headers, std::ostream & ofile )
{
#if defined(FEELPP_HAS_LIBCURL)
    curl_global_init(CURL_GLOBAL_ALL);
    CURL *curl_handle;
    curl_handle = curl_easy_init();
    curl_easy_setopt(curl_handle, CURLOPT_URL, url.c_str() );
    curl_easy_setopt(curl_handle, CURLOPT_HTTPGET, 1L);

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
    if ( auto it = pt.get_optional<std::string>("token") )
        M_token = *it;
    if ( auto it = pt.get_child_optional("file") )
    {
        for( auto const& item : pt.get_child("file") )
            M_fileIds.insert( item.second.get_value<std::string>() );
        if ( M_fileIds.empty() )
            M_fileIds.insert( pt.get<std::string>("file") );
    }

    if ( M_url.empty() )
        M_url = "https://girder.math.unistra.fr";

#if 0
    std::cout << "url: " << M_url << "\n";
    for ( std::string const& fileId : M_fileIds )
        std::cout << "file id: " << fileId << "\n";
#endif

}

bool
RemoteData::Girder::isInit() const
{
    return !M_url.empty();
}
bool
RemoteData::Girder::canDownload() const
{
    return this->isInit() && !M_fileIds.empty();
}


std::vector<std::string>
RemoteData::Girder::download( std::string const& dir ) const
{
    std::vector<std::string> downloadedFileOrFolder;
    if ( M_worldComm.isMasterRank() )
    {
        for ( std::string const& fileId : M_fileIds )
        {
            std::string file = this->downloadFile( fileId, dir );
            if ( !file.empty() )
                downloadedFileOrFolder.push_back( file );
        }
    }
    mpi::broadcast( M_worldComm.globalComm(), downloadedFileOrFolder, M_worldComm.masterRank() );
    return downloadedFileOrFolder;
}

std::string
RemoteData::Girder::downloadFile( std::string const& fileId, std::string const& dir ) const
{
    std::string downloadedFile;
    // get metadata info
    std::string urlFileInfo = M_url+"/api/v1/file/" + fileId;
    std::vector<std::string> headersFileInfo;
    headersFileInfo.push_back("Accept: application/json");
    if ( !M_token.empty() )
        headersFileInfo.push_back( "Girder-Token: "+M_token );
    std::ostringstream omemfile;
    requestHTTPGET( urlFileInfo, headersFileInfo, omemfile );
    // convert to property tree
    std::istringstream istr( omemfile.str() );
    pt::ptree pt;
    pt::read_json(istr, pt);

    // extract info of ptree
    auto itFileName = pt.get_optional<std::string>("name");
    CHECK( itFileName ) << "invalid id : not a file or not exists";
    std::string filename = *itFileName;
    auto itMimeType = pt.get_optional<std::string>("mimeType");

    if ( !fs::exists( dir ) )
        fs::create_directories( dir );

    // download the file
    std::string urlFileDownload = M_url+"/api/v1/file/" + fileId + "/download";
    std::vector<std::string> headersFileDownload;
    if ( itMimeType )
        headersFileDownload.push_back( "Accept: " + *itMimeType );
    if ( !M_token.empty() )
        headersFileDownload.push_back( "Girder-Token: "+M_token );
    std::string filepath = (fs::path(dir)/filename).string();
    std::ofstream ofile( filepath, std::ios::out|std::ios::binary);
    requestHTTPGET( urlFileDownload,headersFileDownload,ofile );
    ofile.close();
    // save metadata
    std::string metadatapath = (fs::path(dir)/(filename+".metadata.json")).string();
    std::ofstream ofileMetadata( metadatapath, std::ios::out);
    ofileMetadata << omemfile.str();
    ofileMetadata.close();

    downloadedFile = filepath;
    return downloadedFile;
}

} // namespace Feel
