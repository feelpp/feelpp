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
#ifndef FEELPP_CORE_REMOTEDATA_HPP
#define FEELPP_CORE_REMOTEDATA_HPP 1

#include <boost/property_tree/json_parser.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>

namespace Feel
{

/**
 * \brief class which manage downloads from an url and download/upload with the Girder database by using the libcurl
 */
struct RemoteData
{
    RemoteData( std::string const& desc, WorldComm const& worldComm = Environment::worldComm() );
    RemoteData( RemoteData const& ) = default;
    RemoteData( RemoteData && ) = default;

    //! return true if enough information are available for download a file
    bool canDownload() const;

    //! download the file
    //! @param dir : the directory where the file is downloaded
    //! @param filename : the filename of the downloaded file
    //! @return : the path of the downloaded file
    std::string download( std::string const& dir = Environment::downloadsRepository(), std::string const& filename = "" ) const;

    class URL
    {
    public :
        URL( std::string const& url, WorldComm const& worldComm = Environment::worldComm() );
        URL( URL const& ) = default;
        URL( URL && ) = default;

        //! return true if the url is valid
        bool isValid() const;

        //! download a file from the url
        //! @param dir : the directory where the file is downloaded
        //! @param filename : the filename of the downloaded file
        //! @return : the path of the downloaded file
        std::string download( std::string const& dir = Environment::downloadsRepository(), std::string const& filename = "" ) const;
    private :
        WorldComm const& M_worldComm;
        std::string M_url;
        std::string M_protocol, M_domain, M_port, M_path, M_query;
    };

    class Github
    {
    public :
        Github( std::string const& desc, WorldComm const& worldComm = Environment::worldComm() );
        Github( Github const& ) = default;
        Github( Github && ) = default;
        void download( std::string const& dir = Environment::downloadsRepository() );
    private :
        WorldComm const& M_worldComm;
        std::string M_owner, M_repo, M_branch, M_path;
    };

    class Girder
    {
    public :
        void download( size_type id, std::string const& dir );
    private :
    };


private :
    boost::optional<URL> M_url;
};

}

#endif
