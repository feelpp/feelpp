//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
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
//! @author Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
//! @date 2024-12-01
//! @copyright 2024 Feel++ Consortium
//!
#include <boost/algorithm/string.hpp>
#include <cpr/cpr.h>
#include <feel/feelcore/remotedata.hpp>
#include <fstream>
#include <regex>


namespace Feel
{

RemoteData::CKAN::CKAN( std::string const& desc, WorldComm& worldComm )
    : M_worldComm( worldComm.shared_from_this() )
{
    std::regex ex("([ ]*)ckan([ ]*):([ ]*)([{])([^]*)([}])");
    std::cmatch what;
    if( !regex_match(desc.c_str(), what, ex) )
        return;

    CHECK( what.size() == 7 ) << "invalid CKAN description format";

    auto resConvertion = convertDescToJson(std::string(what[5].first, what[5].second));
    if (!resConvertion.first)
    {
        if (M_worldComm->isMasterRank())
            std::cout << "CKAN desc has a syntax error: " << desc << "\n";
        return;
    }
    nl::json jsonObj = resConvertion.second;

    if (jsonObj.contains("url"))
        M_url = jsonObj["url"].get<std::string>() + "/api/3/action";
    if (jsonObj.contains("api_key"))
        M_apiKey = jsonObj["api_key"].get<std::string>();
    if (jsonObj.contains("organization"))
        M_organization = jsonObj["organization"].get<std::string>();
    if (jsonObj.contains("dataset"))
        M_dataset = jsonObj["dataset"].get<std::string>();
    if (jsonObj.contains("resource"))
        M_resourceId = jsonObj["resource"].get<std::string>();

    if (M_url.empty())
        M_url = "https://ckan.default.url";

    // Use environment variables for the API key if not specified
    if (M_apiKey.empty())
    {
        // First try FEELPP_CKAN_API_KEY, then fallback to CKAN_API_KEY
        char* env = getenv("FEELPP_CKAN_API_KEY");
        if (env != nullptr && env[0] != '\0')
        {
            M_apiKey = env;
        }
        else
        {
            env = getenv("CKAN_API_KEY");
            if (env != nullptr && env[0] != '\0')
            {
                M_apiKey = env;
            }
        }
    }

    // Use environment variable for organization if not specified
    if (M_organization.empty())
    {
        // First try FEELPP_CKAN_ORGANIZATION, then fallback to CKAN_ORGANIZATION
        char* env = getenv("FEELPP_CKAN_ORGANIZATION");
        if (env != nullptr && env[0] != '\0')
        {
            M_organization = env;
        }
        else
        {
            env = getenv("CKAN_ORGANIZATION");
            if (env != nullptr && env[0] != '\0')
            {
                M_organization = env;
            }
        }
    }
    if ( Environment::isMasterRank() )
    {
        std::cout << fmt::format( "CKAN: url: {}, organization: {}, dataset: {}, resource: {}\n",
                                  M_url, M_organization, M_dataset, M_resourceId ) << "\n";
    }
}

bool RemoteData::CKAN::isInit() const
{
    // Basic initialization check
    return !M_url.empty();
}

bool RemoteData::CKAN::canPerformDatasetOperations() const
{
    return isInit() && !M_apiKey.empty() && !M_organization.empty();
}

bool RemoteData::CKAN::canAccessResource() const
{
    return isInit() && (!M_dataset.empty() || !M_resourceId.empty());
}

// In RemoteData::CKAN class (RemoteData.hpp or RemoteData.cpp)

bool RemoteData::CKAN::canDownload() const
{
    return canAccessResource();
}

bool RemoteData::CKAN::canUpload() const
{
    return canPerformDatasetOperations();
}

std::vector<std::string>
RemoteData::CKAN::download(const std::string& dir) const
{
    std::vector<std::string> downloadedFiles;
    if (!isInit())
        throw std::runtime_error("CKAN remote data object is not initialized.");

    if (M_worldComm->isMasterRank())
    {
        // If we have a specific resource ID, download that resource
        if (!M_resourceId.empty())
        {
            std::string resourceUrl = fmt::format("{}/resource_show", M_url);
            
            // Query resource metadata using CPR
            cpr::Response res = cpr::Get(
                cpr::Url{resourceUrl},
                cpr::Parameters{{"id", M_resourceId}},
                cpr::Header{{"Authorization", M_apiKey}},
                cpr::VerifySsl{false}
            );

            if (res.status_code != 200)
            {
                throw std::runtime_error(fmt::format("CKAN resource query failed: HTTP {}, {}", res.status_code, res.text));
            }

            // Parse JSON response
            nl::json jsonResponse = nl::json::parse(res.text);
            if (!jsonResponse.contains("result"))
                throw std::runtime_error("Invalid CKAN resource response structure.");

            std::string downloadUrl = jsonResponse["result"].value("url", "");
            std::string filename = jsonResponse["result"].value("name", "downloaded_resource");
            if (downloadUrl.empty())
                throw std::runtime_error("CKAN resource URL is empty.");

            fs::path filepath = fs::path(dir) / filename;
            
            // Download the actual file content
            cpr::Response downloadRes = cpr::Get(
                cpr::Url{downloadUrl},
                cpr::Header{{"Authorization", M_apiKey}},
                cpr::VerifySsl{false}
            );

            if (downloadRes.status_code != 200)
            {
                throw std::runtime_error(fmt::format("CKAN resource download failed: HTTP {}, {}", downloadRes.status_code, downloadRes.text));
            }

            // Write the file content to disk
            std::ofstream file(filepath.string(), std::ios::binary);
            file.write(downloadRes.text.data(), downloadRes.text.size());
            file.close();

            if (!file.good())
            {
                throw std::runtime_error(fmt::format("Failed to write downloaded file: {}", filepath.string()));
            }

            downloadedFiles.push_back(filepath.string());
        }
        // If we have a dataset, download all resources from that dataset
        else if (!M_dataset.empty())
        {
            std::string packageUrl = fmt::format("{}/package_show", M_url);
            
            // Query dataset metadata using CPR
            cpr::Response res = cpr::Get(
                cpr::Url{packageUrl},
                cpr::Parameters{{"id", M_dataset}},
                cpr::Header{{"Authorization", M_apiKey}},
                cpr::VerifySsl{false}
            );

            if (res.status_code != 200)
            {
                throw std::runtime_error(fmt::format("CKAN dataset query failed: HTTP {}, {}", res.status_code, res.text));
            }

            // Parse JSON response
            nl::json jsonResponse = nl::json::parse(res.text);
            if (!jsonResponse.contains("result") || !jsonResponse["result"].contains("resources"))
                throw std::runtime_error("Invalid CKAN dataset response structure.");

            // Download each resource in the dataset
            auto resources = jsonResponse["result"]["resources"];
            for (const auto& resource : resources)
            {
                std::string downloadUrl = resource.value("url", "");
                std::string filename = resource.value("name", "downloaded_resource");
                
                if (downloadUrl.empty())
                {
                    std::cerr << "Warning: Skipping resource with empty URL: " << filename << std::endl;
                    continue;
                }

                fs::path filepath = fs::path(dir) / filename;
                
                // Download the actual file content
                cpr::Response downloadRes = cpr::Get(
                    cpr::Url{downloadUrl},
                    cpr::Header{{"Authorization", M_apiKey}},
                    cpr::VerifySsl{false}
                );

                if (downloadRes.status_code != 200)
                {
                    std::cerr << "Warning: Failed to download resource " << filename 
                              << ": HTTP " << downloadRes.status_code << std::endl;
                    continue;
                }

                // Write the file content to disk
                std::ofstream file(filepath.string(), std::ios::binary);
                file.write(downloadRes.text.data(), downloadRes.text.size());
                file.close();

                if (!file.good())
                {
                    std::cerr << "Warning: Failed to write file: " << filepath.string() << std::endl;
                    continue;
                }

                downloadedFiles.push_back(filepath.string());
                std::cout << "Downloaded: " << filename << " (" << downloadRes.text.size() << " bytes)" << std::endl;
            }
        }
        else
        {
            throw std::runtime_error("No resource ID or dataset specified for CKAN download.");
        }
    }

    M_worldComm->barrier();
    return downloadedFiles;
}

std::vector<std::string>
RemoteData::CKAN::upload(const std::string& dataPath, const std::string& datasetId, bool sync) const
{
    CHECK(isInit()) << "CKAN remote data object is not initialized.";
    CHECK(fs::exists(dataPath)) << fmt::format("Data path '{}' does not exist.", dataPath);

    std::vector<std::string> uploadedResources;
    if (M_worldComm->isMasterRank())
    {
        std::string url = fmt::format("{}/resource_create", M_url);
        cpr::Header headers = {
            {"Authorization", M_apiKey}//,
            //{"Accept", "application/json"}
        };

        // Prepare file payload
        fs::path filepath(dataPath);
        std::ifstream fileStream(dataPath, std::ios::binary | std::ios::ate);
        std::streamsize fileSize = fileStream.tellg();
        fileStream.seekg(0, std::ios::beg);
        std::vector<char> fileData(fileSize);
        fileStream.read(fileData.data(), fileSize);

        // HTTP POST request
        cpr::Multipart multipart = {
            {"package_id", datasetId},
            {"name", filepath.filename().string()},
            {"upload", cpr::Buffer(fileData.begin(), fileData.end(), filepath.filename().string())}
        };

        cpr::Response res = cpr::Post(
            cpr::Url{url},
            headers,
            multipart,
            cpr::VerifySsl{false}
        );
        std::cout << fmt::format("CKAN upload response: {}", res.text) << "\n";
        if (res.status_code != 200)
        {
            nl::json errorResponse = nl::json::parse(res.text);
            std::string errorMsg = errorResponse.contains("error") ? errorResponse["error"].dump() : "Unknown error";
            throw std::runtime_error(fmt::format("CKAN upload failed: {}", errorMsg));
        }

        // Parse the JSON response for the uploaded resource details
        nl::json jsonResponse = nl::json::parse(res.text);
        uploadedResources.push_back(jsonResponse["result"].value("id", ""));

        std::cout << fmt::format("Uploaded resource with ID: {}", uploadedResources[0]) << "\n";
    }

    if (sync)
        mpi::broadcast(M_worldComm->globalComm(), uploadedResources, M_worldComm->masterRank());

    return uploadedResources;
}
std::string
RemoteData::CKAN::createDataset(const std::string& name, const std::string& organization, const std::string& description) const
{
    if (!canPerformDatasetOperations())
    {
        throw std::runtime_error("CKAN is not properly initialized for dataset operations.");
    }

    std::string url = fmt::format("{}/package_create", M_url);
    cpr::Header headers = {
        {"Authorization", M_apiKey},
        {"Content-Type", "application/json"}
    };

    // Prepare the dataset metadata
    nl::json datasetMetadata = {
        {"name", name},
        {"owner_org", organization},
        {"title", name},
        {"notes", description}
    };

    // HTTP POST request
    cpr::Response res = cpr::Post(
        cpr::Url{url},
        headers,
        cpr::Payload{
            {"name", name},
            {"owner_org", organization},
            {"title", name},
            {"notes", description}
        },
        cpr::VerifySsl{false}
    );

    if (res.status_code != 200)
    {
        nl::json errorResponse = nl::json::parse(res.text);
        std::string errorMsg = errorResponse.contains("error") ? errorResponse["error"].dump() : "Unknown error";
        throw std::runtime_error(fmt::format("CKAN dataset creation failed: {}", errorMsg));
    }

    // Parse the JSON response for the dataset ID
    nl::json jsonResponse = nl::json::parse(res.text);
    return jsonResponse["result"].value("id", "");
}

// In RemoteData::CKAN class

nl::json RemoteData::CKAN::createDataset(const std::string& datasetName) const
{
    // Construct the URL for creating a dataset
    std::string url = fmt::format("{}/package_create", M_url);

    // Prepare headers
    cpr::Header headers = {
        {"Authorization", M_apiKey},
        {"Content-Type", "application/json"}
    };
    nl::json payload = {
        {"name", datasetName},
        {"notes", "Created by Cemosis"},
        {"owner_org", M_organization}
    };
    // Send the POST request
    cpr::Response res = cpr::Post(
        cpr::Url{url},
        headers,
        cpr::Body{payload.dump()},
        cpr::VerifySsl{false}
    );

    if (res.status_code != 200)
    {
        throw std::runtime_error(fmt::format("Failed to create dataset on CKAN: {}, {}, res: {}", M_url, M_organization, res.text));
    }
    // Parse and return the JSON response
    nl::json jsonResponse = nl::json::parse(res.text);
    return jsonResponse["result"];
}

bool RemoteData::CKAN::deleteDataset(const std::string& datasetId) const
{
    if (!canPerformDatasetOperations())
    {
        throw std::runtime_error("CKAN is not properly initialized for dataset operations.");
    }
    // Construct the URL for deleting a dataset
    std::string url = fmt::format("{}/dataset_purge", M_url);

    // Prepare headers
    cpr::Header headers = {
        {"Authorization", M_apiKey},
        {"Content-Type", "application/json"}
    };

    // Send the POST request
    cpr::Response res = cpr::Post(
        cpr::Url{url},
        headers,
        cpr::Payload{
            {"id", datasetId}
        },
        cpr::VerifySsl{false}
    );

    return res.status_code == 200;
}

std::tuple<std::vector<std::shared_ptr<RemoteData::FolderInfo>>, std::vector<std::shared_ptr<RemoteData::ItemInfo>>, std::vector<std::shared_ptr<RemoteData::FileInfo>>>
RemoteData::CKAN::contents() const
{
    using FolderInfo = RemoteData::FolderInfo;
    using ItemInfo = RemoteData::ItemInfo;
    using FileInfo = RemoteData::FileInfo;
    
    std::vector<std::shared_ptr<FolderInfo>> folders;
    std::vector<std::shared_ptr<ItemInfo>> items;
    std::vector<std::shared_ptr<FileInfo>> files;

    if ( !canDownload() )
        return std::make_tuple( folders, items, files );

    try 
    {
        // Get dataset show API endpoint
        std::string url = M_url + "/package_show";
        
        // Prepare POST data
        nl::json requestData;
        requestData["id"] = M_dataset;
        
        cpr::Response r = cpr::Post(
            cpr::Url{url},
            cpr::Body{requestData.dump()},
            cpr::Header{{"Authorization", M_apiKey}, {"Content-Type", "application/json"}}
        );

        if ( r.status_code == 200 )
        {
            auto response = nl::json::parse( r.text );
            if ( response["success"].get<bool>() )
            {
                auto dataset = response["result"];
                
                // Create a folder representing the dataset itself
                auto datasetFolder = std::make_shared<FolderInfo>( 
                    dataset.value("title", dataset.value("name", M_dataset)), 
                    dataset.value("id", M_dataset), 
                    dataset.value("num_resources", 0)
                );
                folders.push_back( datasetFolder );

                // Process resources as files
                if ( dataset.contains("resources") && dataset["resources"].is_array() )
                {
                    for ( const auto& resource : dataset["resources"] )
                    {
                        auto fileInfo = std::make_shared<FileInfo>(
                            resource.value("name", resource.value("id", "unknown")),
                            resource.value("id", ""),
                            resource.value("size", 0)
                        );
                        files.push_back( fileInfo );
                    }
                }
            }
        }
    }
    catch ( const std::exception& e )
    {
        Feel::cout << "Error getting CKAN dataset contents: " << e.what() << std::endl;
    }

    return std::make_tuple( folders, items, files );
}

std::vector<std::string>
RemoteData::CKAN::upload(const std::string& dataPath, const std::string& parentId) const
{
    // Use the dataset as the parent ID context for CKAN uploads
    return upload(dataPath, parentId.empty() ? M_dataset : parentId, true);
}

} // namespace Feel