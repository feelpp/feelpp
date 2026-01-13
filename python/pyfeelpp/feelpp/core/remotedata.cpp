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
#include <feel/feelpython/pybind11/functional.h>

#include<feel/feelcore/environment.hpp>
#include<feel/feelcore/remotedata.hpp>
#include<feel/feelcore/json.hpp>

namespace py = pybind11;
void
bindRemoteData( py::module & m )
{
    using namespace Feel;
    
    // Bind StatusRequestHTTP class
    py::class_<StatusRequestHTTP>(m, "StatusRequestHTTP")
        .def(py::init<bool, std::string const&>(),
             py::arg("s"), py::arg("m")=std::string(""),
             "Initialize StatusRequestHTTP")
        .def(py::init<bool, uint16_type, std::string const&>(),
             py::arg("s"), py::arg("c"), py::arg("m")=std::string(""),
             "Initialize StatusRequestHTTP with code")
        .def("success", &StatusRequestHTTP::success, "check if request was successful")
        .def("code", &StatusRequestHTTP::code, "get HTTP status code")
        .def("msg", &StatusRequestHTTP::msg, "get status message")
        ;
    
    // Bind utility functions
    m.def("requestHTTPGET", &requestHTTPGET,
          py::arg("url"), py::arg("headers"), py::arg("ofile"),
          py::arg("timeout")=5000, py::arg("max_retries")=3, py::arg("backoff_delay")=1000,
          "make HTTP GET request");
    
    m.def("requestHTTPPOST", static_cast<StatusRequestHTTP(*)(const std::string&, const std::vector<std::string>&, std::ostream&, int, int, int)>(&requestHTTPPOST),
          py::arg("url"), py::arg("headers"), py::arg("ofile"),
          py::arg("timeout")=5000, py::arg("max_retries")=3, py::arg("backoff_delay")=1000,
          "make HTTP POST request");
    
    m.def("requestDownloadURL", &requestDownloadURL,
          py::arg("url"), py::arg("ofile"),
          py::arg("timeout")=5000, py::arg("max_retries")=3, py::arg("backoff_delay")=1000,
          "download from URL");
    
    m.def("convertDescToJson", &convertDescToJson,
          py::arg("desc"),
          "convert description to JSON");
    
    // Bind RemoteDataProgress class
    py::enum_<RemoteDataProgress::Operation>(m, "RemoteDataProgressOperation")
        .value("UPLOAD", RemoteDataProgress::Operation::UPLOAD)
        .value("DOWNLOAD", RemoteDataProgress::Operation::DOWNLOAD)
        ;
    
    py::enum_<RemoteDataProgress::Level>(m, "RemoteDataProgressLevel")
        .value("QUIET", RemoteDataProgress::Level::QUIET)
        .value("NORMAL", RemoteDataProgress::Level::NORMAL)
        .value("VERBOSE", RemoteDataProgress::Level::VERBOSE)
        .value("DEBUG", RemoteDataProgress::Level::DEBUG)
        ;
    
    py::class_<RemoteDataProgress>(m, "RemoteDataProgress")
        .def(py::init<RemoteDataProgress::Operation, RemoteDataProgress::Level>(),
             py::arg("op"), py::arg("level") = RemoteDataProgress::Level::NORMAL,
             "Initialize RemoteDataProgress")
        .def("setLevel", &RemoteDataProgress::setLevel, "set verbosity level")
        .def("isQuiet", &RemoteDataProgress::isQuiet, "check if quiet mode")
        .def("isNormal", &RemoteDataProgress::isNormal, "check if normal verbosity")
        .def("isVerbose", &RemoteDataProgress::isVerbose, "check if verbose mode")
        .def("isDebug", &RemoteDataProgress::isDebug, "check if debug mode")
        .def("formatSize", &RemoteDataProgress::formatSize, "format file size in human readable format")
        ;
    
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
        .def(py::init<std::string const&,worldcomm_ptr_t const&>(),py::arg("desc"),py::arg("worldComm"),"Initialize the RemoteData handler")
        .def(py::init([](std::string const& desc) {
            return RemoteData(desc, Environment::worldCommPtr());
        }), py::arg("desc"), "Initialize the RemoteData handler with default worldComm")
        .def("worldComm", &RemoteData::worldComm, "get the worldComm" )
        .def("canDownload", &RemoteData::canDownload, "returns true if data/ressource can be downloaded, false otherwise" )
        .def("canUpload", &RemoteData::canUpload, "returns true if data/ressource can be uploaded, false otherwise" )
        .def("download", [](RemoteData const& self, std::string dir, std::string filename) {
                 if (dir.empty()) dir = Environment::downloadsRepository();
                 return self.download(dir, filename);
             }, py::arg("dir")="", py::arg("filename")="",
             "download the requested data/ressource" )
        .def("download", static_cast<std::vector<std::string> (RemoteData::*)( std::string const&, std::string const&, int) const>(&RemoteData::download),
             py::arg("dir"), py::arg("filename"), py::arg("timeout"),
             "download the requested data/ressource with timeout" )
        .def("upload", static_cast<std::vector<std::string> (RemoteData::*)( std::string const&, std::string const&, bool) const>(&RemoteData::upload),
             py::arg("path"),
             py::arg("parentId")=std::string(""),
             py::arg("sync")=true,
             "upload the requested data/ressource" )
        .def("upload", static_cast<std::vector<std::string> (RemoteData::*)( std::string const&, std::string const&, bool, int) const>(&RemoteData::upload),
             py::arg("path"), py::arg("parentId"), py::arg("sync"), py::arg("timeout"),
             "upload the requested data/ressource with timeout" )
        .def("upload", static_cast<std::vector<std::vector<std::string>> (RemoteData::*)( std::vector<std::pair<std::string, std::string>> const&, bool) const>(&RemoteData::upload),
             py::arg("dataToUpload"), py::arg("sync")=true,
             "upload multiple data/ressources" )
        .def("replaceFile", static_cast<void (RemoteData::*)( std::string const&, std::string const&) const>(&RemoteData::replaceFile),
             py::arg("filePath"), py::arg("fileId"),
             "replace contents of a file" )
        .def("replaceFile", static_cast<void (RemoteData::*)( std::vector<std::pair<std::string, std::string>> const&) const>(&RemoteData::replaceFile),
             py::arg("filesToReplace"),
             "replace contents of multiple files" )
        .def("createFolder", &RemoteData::createFolder,
             py::arg("folderPath"), py::arg("parentId")=std::string(""), py::arg("sync")=true,
             "create folders hierarchy on remote storage" )
        .def("createItem", &RemoteData::createItem,
             py::arg("itemPath"), py::arg("parentId"), py::arg("sync")=true,
             "create item on remote storage" )
        .def("createDataset", &RemoteData::createDataset,
             py::arg("datasetName"),
             "create dataset" )
        .def("deleteDataset", &RemoteData::deleteDataset,
             py::arg("datasetId"),
             "delete dataset" )
        .def("resourceLookup", &RemoteData::resourceLookup,
             py::arg("path"), py::arg("token")=std::string(""),
             "lookup a resource by path" )
        .def("deleteResource", &RemoteData::deleteResource,
             py::arg("resourceId"), py::arg("token")=std::string(""),
             "delete a resource by id" )
        .def("listOrganizations", &RemoteData::listOrganizations,
             "list organizations available on the remote data platform" )
        .def( "contents", static_cast<RemoteData::ContentsInfo (RemoteData::*)() const>(&RemoteData::contents), 
              "get the data/ressource contents information" )
        .def( "contents", static_cast<RemoteData::ContentsInfo (RemoteData::*)(RemoteDataProgress&) const>(&RemoteData::contents), 
              py::arg("progress"), "get the data/ressource contents information with progress control" )
        ;
    
    // Bind RemoteData::URL class
    py::class_<RemoteData::URL>(m, "RemoteDataURL")
        .def(py::init([](std::string const& url, py::object worldComm) {
                 if (worldComm.is_none())
                     return new RemoteData::URL(url, Environment::worldComm());
                 return new RemoteData::URL(url, worldComm.cast<WorldComm&>());
             }), py::arg("url"), py::arg("worldComm")=py::none(),
             "Initialize URL handler")
        .def("isValid", &RemoteData::URL::isValid, "return true if the URL is valid")
        .def("download", [](RemoteData::URL const& self, std::string dir, std::string filename) {
                 if (dir.empty()) dir = Environment::downloadsRepository();
                 return self.download(dir, filename);
             }, py::arg("dir")="", py::arg("filename")="",
             "download a file from the URL")
        ;
    
    // Bind RemoteData::Github class
    py::class_<RemoteData::Github>(m, "RemoteDataGithub")
        .def(py::init([](std::string const& desc, py::object worldComm) {
                 if (worldComm.is_none())
                     return new RemoteData::Github(desc, Environment::worldComm());
                 return new RemoteData::Github(desc, worldComm.cast<WorldComm&>());
             }), py::arg("desc"), py::arg("worldComm")=py::none(),
             "Initialize Github handler")
        .def("isInit", &RemoteData::Github::isInit, "return true if Github is initialized")
        .def("canDownload", &RemoteData::Github::canDownload, "return true if can download")
        .def("canUpload", &RemoteData::Github::canUpload, "return true if can upload")
        .def("download", [](RemoteData::Github const& self, std::string dir) {
                 if (dir.empty()) dir = Environment::downloadsRepository();
                 return self.download(dir);
             }, py::arg("dir")="",
             "download file/folder from Github")
        ;
    
    // Bind RemoteData::Girder class
    py::class_<RemoteData::Girder>(m, "RemoteDataGirder")
        .def(py::init([](std::string const& desc, py::object worldComm) {
                 if (worldComm.is_none())
                     return new RemoteData::Girder(desc, Environment::worldComm());
                 return new RemoteData::Girder(desc, worldComm.cast<WorldComm&>());
             }), py::arg("desc"), py::arg("worldComm")=py::none(),
             "Initialize Girder handler")
        .def("setFolderIds", &RemoteData::Girder::setFolderIds, "set folder ids to only one folder id")
        .def("isInit", &RemoteData::Girder::isInit, "return true if Girder is initialized")
        .def("canDownload", &RemoteData::Girder::canDownload, "return true if can download")
        .def("canUpload", &RemoteData::Girder::canUpload, "return true if can upload")
        .def("download", [](RemoteData::Girder const& self, std::string dir) {
                 if (dir.empty()) dir = Environment::downloadsRepository();
                 return self.download(dir);
             }, py::arg("dir")="",
             "download file/folder from Girder")
        .def("download", static_cast<std::vector<std::string> (RemoteData::Girder::*)( std::string const&, int) const>(&RemoteData::Girder::download),
             py::arg("dir"), py::arg("timeout"),
             "download file/folder from Girder with timeout")
        .def("download", static_cast<std::vector<std::string> (RemoteData::Girder::*)( std::string const&, std::string const&) const>(&RemoteData::Girder::download),
             py::arg("dir"), py::arg("path"),
             "download file/folder/item from Girder by path")
        // TODO: Fix Girder download with timeout and path - method declared but not implemented
        // .def("download", static_cast<std::vector<std::string> (RemoteData::Girder::*)( std::string const&, std::string const&, int) const>(&RemoteData::Girder::download),
        //      py::arg("dir"), py::arg("path"), py::arg("timeout"),
        //      "download file/folder/item from Girder by path with timeout")
        .def("upload", static_cast<std::vector<std::string> (RemoteData::Girder::*)( std::string const&, std::string const&, bool) const>(&RemoteData::Girder::upload),
             py::arg("dataPath"), py::arg("parentId")=std::string(""), py::arg("sync")=true,
             "upload data on Girder")
        .def("upload", static_cast<std::vector<std::string> (RemoteData::Girder::*)( std::string const&, std::string const&, bool, int) const>(&RemoteData::Girder::upload),
             py::arg("dataPath"), py::arg("parentId"), py::arg("sync"), py::arg("timeout"),
             "upload data on Girder with timeout")
        .def("upload", static_cast<std::vector<std::vector<std::string>> (RemoteData::Girder::*)( std::vector<std::pair<std::string, std::string>> const&, bool) const>(&RemoteData::Girder::upload),
             py::arg("dataToUpload"), py::arg("sync")=true,
             "upload multiple data on Girder")
        .def("replaceFile", static_cast<void (RemoteData::Girder::*)( std::string const&, std::string const&) const>(&RemoteData::Girder::replaceFile),
             py::arg("filePath"), py::arg("fileId"),
             "replace contents of a file")
        .def("replaceFile", static_cast<void (RemoteData::Girder::*)( std::vector<std::pair<std::string, std::string>> const&) const>(&RemoteData::Girder::replaceFile),
             py::arg("filesToReplace"),
             "replace contents of files")
        .def("createFolder", &RemoteData::Girder::createFolder,
             py::arg("folderPath"), py::arg("parentId")=std::string(""), py::arg("sync")=true,
             "create folders hierarchy on Girder")
        .def("createItem", &RemoteData::Girder::createItem,
             py::arg("itemPath"), py::arg("parentId"), py::arg("sync")=true,
             "create item on Girder")
        .def("contents", static_cast<std::tuple<std::vector<std::shared_ptr<RemoteData::FolderInfo>>, std::vector<std::shared_ptr<RemoteData::ItemInfo>>, std::vector<std::shared_ptr<RemoteData::FileInfo>>> (RemoteData::Girder::*)() const>(&RemoteData::Girder::contents), 
             "get contents of remote data")
        .def("contents", static_cast<std::tuple<std::vector<std::shared_ptr<RemoteData::FolderInfo>>, std::vector<std::shared_ptr<RemoteData::ItemInfo>>, std::vector<std::shared_ptr<RemoteData::FileInfo>>> (RemoteData::Girder::*)(RemoteDataProgress&) const>(&RemoteData::Girder::contents), 
             py::arg("progress"), "get contents of remote data with progress control")
        .def("resourceLookup", static_cast<nl::json (RemoteData::Girder::*)(const std::string&, const std::string&) const>(&RemoteData::Girder::resourceLookup),
             py::arg("path"), py::arg("token")=std::string(""),
             "lookup a resource by path")
        .def("resourceLookup", static_cast<nl::json (RemoteData::Girder::*)(const std::string&, const std::string&, const RemoteDataProgress&) const>(&RemoteData::Girder::resourceLookup),
             py::arg("path"), py::arg("token"), py::arg("progress"),
             "lookup a resource by path with progress control")
        .def("deleteResource", &RemoteData::Girder::deleteResource,
             py::arg("resourceId"), py::arg("token")=std::string(""),
             "delete a resource")
        ;
    
    // Bind RemoteData::CKAN class
    py::class_<RemoteData::CKAN>(m, "RemoteDataCKAN")
        .def("isInit", &RemoteData::CKAN::isInit, "return true if CKAN is initialized")
        .def("canDownload", &RemoteData::CKAN::canDownload, "return true if can download")
        .def("canUpload", &RemoteData::CKAN::canUpload, "return true if can upload")
        .def("download", [](RemoteData::CKAN const& self, std::string dir) {
                 if (dir.empty()) dir = Environment::downloadsRepository();
                 return self.download(dir);
             }, py::arg("dir")="",
             "download from CKAN")
        .def("download", static_cast<std::vector<std::string> (RemoteData::CKAN::*)( std::string const&, int) const>(&RemoteData::CKAN::download),
             py::arg("dir"), py::arg("timeout"),
             "download from CKAN with timeout")
        .def("contents", &RemoteData::CKAN::contents, "get contents of CKAN data")
        .def("upload", static_cast<std::vector<std::string> (RemoteData::CKAN::*)( std::string const&, std::string const&) const>(&RemoteData::CKAN::upload),
             py::arg("dataPath"), py::arg("parentId")=std::string(""),
             "upload data to CKAN")
        .def("upload", static_cast<std::vector<std::string> (RemoteData::CKAN::*)( std::string const&, std::string const&, int) const>(&RemoteData::CKAN::upload),
             py::arg("dataPath"), py::arg("parentId"), py::arg("timeout"),
             "upload data to CKAN with timeout")
        .def("upload", static_cast<std::vector<std::string> (RemoteData::CKAN::*)( std::string const&, std::string const&, bool) const>(&RemoteData::CKAN::upload),
             py::arg("dataPath"), py::arg("datasetId"), py::arg("sync"),
             "upload data to CKAN dataset")
        // TODO: Fix CKAN replaceResource binding - method not implemented
        // .def("replaceResource", &RemoteData::CKAN::replaceResource,
        //      py::arg("resourcePath"), py::arg("resourceId"),
        //      "replace a resource")
        // TODO: Fix CKAN createResource binding - method not implemented  
        // .def("createResource", &RemoteData::CKAN::createResource,
        //      py::arg("name"), py::arg("description"), py::arg("parentId"),
        //      "create a resource")
        .def("createDataset", static_cast<std::string (RemoteData::CKAN::*)( std::string const&, std::string const&, std::string const&) const>(&RemoteData::CKAN::createDataset),
             py::arg("name"), py::arg("organization"), py::arg("description"),
             "create a dataset")
        .def("createDataset", static_cast<nl::json (RemoteData::CKAN::*)( std::string const&) const>(&RemoteData::CKAN::createDataset),
             py::arg("datasetName"),
             "create a dataset by name")
        .def("deleteDataset", &RemoteData::CKAN::deleteDataset,
             py::arg("datasetId"),
             "delete a dataset")
        // TODO: Fix CKAN resourceLookup binding - symbol not found
        // .def("resourceLookup", &RemoteData::CKAN::resourceLookup,
        //      py::arg("pathOrId"),
        //      "lookup a resource")
        // TODO: Fix CKAN deleteResource binding - method not implemented
        // .def("deleteResource", &RemoteData::CKAN::deleteResource,
        //      py::arg("resourceId"),
        //      "delete a resource")
        .def("listOrganizations", &RemoteData::CKAN::listOrganizations,
             "list organizations")
        ;
    
}
