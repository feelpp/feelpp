//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! Copyright (C) 2017-present Feel++ Consortium
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
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 07 Aug 2018
//! @copyright 2018 Feel++ Consortium
//!
#include <feel/feelpython/pybind11/pybind11.h>
#include <feel/feelpython/pybind11/stl.h>


#include<feel/feelcore/environment.hpp>
#include<feel/feelcore/remotedata.hpp>

namespace py = pybind11;
void
bindRemoteData( py::module & m )
{
    using namespace Feel;
    py::class_<RemoteData::ContentsInfo>(m,"ContentsInfo")
        .def( py::init<>() )
        .def( py::init<std::tuple<std::vector<std::shared_ptr<RemoteData::FolderInfo>>,std::vector<std::shared_ptr<RemoteData::ItemInfo>>,std::vector<std::shared_ptr<RemoteData::FileInfo>>>>(), "set ContentsInfo from tuple" )
        .def( "folderInfo", &RemoteData::ContentsInfo::folderInfo, "get remote contents folderInfo" )
        .def( "itemInfo", &RemoteData::ContentsInfo::itemInfo, "get remote contents itemInfo" )
        .def( "fileInfo", &RemoteData::ContentsInfo::fileInfo, "get remote contents fileInfo" )
        ;
    
    py::class_<RemoteData::FolderInfo>(m,"FolderInfo")
        .def( py::init<>() )
        .def( py::init<std::string const&,std::string const&,size_type>(),
              py::arg("name")="",
              py::arg("id")="",
              py::arg("size")=invalid_size_type_value,
              "Initialize FolderInfo" )
        .def( "print", &RemoteData::FolderInfo::print, "print FolderInfo to ostringstream" )
        .def( "name", &RemoteData::FolderInfo::name, "get FolderInfo name" )
        .def( "id", &RemoteData::FolderInfo::id, "get FolderInfo id" )
        .def( "size", &RemoteData::FolderInfo::size, "get FolderInfo size" )
        ;

    py::class_<RemoteData::ItemInfo>(m,"ItemInfo")
        .def( py::init<>() )
        .def( py::init<std::string const&,std::string const&,size_type>(),
              py::arg("name")="",
              py::arg("id")="",
              py::arg("size")=invalid_size_type_value,
              "Initialize FolderInfo" )
        .def( "print", &RemoteData::ItemInfo::print, "print ItemInfo to ostringstream" )
        .def( "name", &RemoteData::ItemInfo::name, "get ItemInfo name" )
        .def( "id", &RemoteData::ItemInfo::id, "get ItemInfo id" )
        .def( "size", &RemoteData::ItemInfo::size, "get ItemInfo size" )
        ;

    py::class_<RemoteData::FileInfo>(m,"FileInfo")
        .def( py::init<>() )
        .def( py::init<std::string const&,std::string const&,size_type>(),
              py::arg("name")="",
              py::arg("id")="",
              py::arg("size")=invalid_size_type_value,
              "Initialize FolderInfo" )
        .def( "print", &RemoteData::FileInfo::print, "print FileInfo to ostringstream" )
        .def( "name", &RemoteData::FileInfo::name, "get FileInfo name" )
        .def( "id", &RemoteData::FileInfo::id, "get FileInfo id" )
        .def( "size", &RemoteData::FileInfo::size, "get FileInfo size" )
        .def( "mimeType", &RemoteData::FileInfo::mimeType, "get FileInfo mimeType" )
        .def( "checksum", &RemoteData::FileInfo::checksum, "get FileInfo checksum" )
        .def( "checksumType", &RemoteData::FileInfo::checksumType, "get FileInfo checksumType" )
        .def( "setChecksum", &RemoteData::FileInfo::setChecksum, "set FileInfo checksum type and value" )
        .def( "setMimeType", &RemoteData::FileInfo::setMimeType, "set FileInfo mimetype" )
        ;
        
    py::class_<RemoteData>(m,"RemoteData")
        .def(py::init<std::string const&,worldcomm_ptr_t const&>(),py::arg("desc"),py::arg("worldComm")=Environment::worldCommPtr(),"Initialize the RemoteData handler")
        .def("worldComm", &RemoteData::worldComm, "get the worldComm" )
        .def("canDownload", &RemoteData::canDownload, "returns true if data/ressource can be downloaded, false otherwise" )
        .def("canUpload", &RemoteData::canUpload, "returns true if data/ressource can be uploaded, false otherwise" )
        .def("download", &RemoteData::download,
             py::arg("dir")=Environment::downloadsRepository(),
             py::arg("filename")=std::string(""),
             "download the requested data/ressource" )
        .def("upload", static_cast<std::vector<std::string> (RemoteData::*)( std::string const&, std::string const&, bool sync) const>(&RemoteData::upload),
             py::arg("path"),
             py::arg("parentId")=std::string(""),
             py::arg("sync")=true,
             "upload the requested data/ressource" )
        .def( "contents", &RemoteData::contents, "get the data/ressource contents information" )
        ;
    
}
