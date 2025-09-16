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
    {
        // Check environment variable for URL if not specified in JSON
        char* env = getenv("FEELPP_CKAN_URL");
        if (env != nullptr && env[0] != '\0')
        {
            M_url = std::string(env) + "/api/3/action";
        }
        else
        {
            env = getenv("CKAN_URL");
            if (env != nullptr && env[0] != '\0')
            {
                M_url = std::string(env) + "/api/3/action";
            }
            else
            {
                M_url = "https://ckan.default.url";
            }
        }
    }

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
        // Determine progress level
        RemoteDataProgress::Level progressLevel = RemoteDataProgress::Level::NORMAL;
        if (Environment::vm().count("quiet")) {
            progressLevel = RemoteDataProgress::Level::QUIET;
        } else if (Environment::vm().count("debug")) {
            progressLevel = RemoteDataProgress::Level::DEBUG;
        } else if (Environment::vm().count("verbose") || Environment::vm().count("progress")) {
            progressLevel = RemoteDataProgress::Level::VERBOSE;
        }
        
        RemoteDataProgress progress(RemoteDataProgress::Operation::DOWNLOAD, progressLevel);
        
        // If we have a specific resource ID, download that resource
        if (!M_resourceId.empty())
        {
            progress.startOperation("Download CKAN resource " + M_resourceId);
            
            std::string resourceUrl = fmt::format("{}/resource_show", M_url);
            
            if (progress.showDebugOutput()) {
                progress.debug("Resource metadata URL: " + resourceUrl);
            }
            
            // Query resource metadata using CPR
            cpr::Response res = cpr::Get(
                cpr::Url{resourceUrl},
                cpr::Parameters{{"id", M_resourceId}},
                cpr::Header{{"Authorization", M_apiKey}},
                cpr::VerifySsl{false},
                cpr::Timeout{30000}  // 30 second default timeout
            );

            if (res.status_code != 200)
            {
                progress.error(fmt::format("CKAN resource query failed: HTTP {}, {}", res.status_code, res.text));
                throw std::runtime_error(fmt::format("CKAN resource query failed: HTTP {}, {}", res.status_code, res.text));
            }

            // Parse JSON response
            nl::json jsonResponse = nl::json::parse(res.text);
            if (progress.showDebugOutput()) {
                progress.debug("Resource metadata response: " + jsonResponse.dump(2));
            }
            
            if (!jsonResponse.contains("result"))
            {
                progress.error("Invalid CKAN resource response structure.");
                throw std::runtime_error("Invalid CKAN resource response structure.");
            }

            std::string downloadUrl = jsonResponse["result"].value("url", "");
            std::string filename = jsonResponse["result"].value("name", "downloaded_resource");
            
            // Handle null or missing size values safely
            std::streamsize fileSize = 0;
            if (jsonResponse["result"].contains("size") && !jsonResponse["result"]["size"].is_null()) {
                fileSize = jsonResponse["result"]["size"].get<std::streamsize>();
            }
            
            if (downloadUrl.empty())
            {
                progress.error("CKAN resource URL is empty.");
                throw std::runtime_error("CKAN resource URL is empty.");
            }

            fs::path filepath = fs::path(dir) / filename;
            
            progress.startFile(filename, fileSize);
            
            if (progress.showDebugOutput()) {
                progress.debug("Download URL: " + downloadUrl);
            }
            
            // Download the actual file content
            cpr::Response downloadRes = cpr::Get(
                cpr::Url{downloadUrl},
                cpr::Header{{"Authorization", M_apiKey}},
                cpr::VerifySsl{false},
                cpr::Timeout{30000}  // 30 second default timeout
            );

            if (downloadRes.status_code != 200)
            {
                progress.error(fmt::format("CKAN resource download failed: HTTP {}, {}", downloadRes.status_code, downloadRes.text));
                throw std::runtime_error(fmt::format("CKAN resource download failed: HTTP {}, {}", downloadRes.status_code, downloadRes.text));
            }

            // Write the file content to disk
            std::ofstream file(filepath.string(), std::ios::binary);
            file.write(downloadRes.text.data(), downloadRes.text.size());
            file.close();

            if (!file.good())
            {
                progress.error(fmt::format("Failed to write downloaded file: {}", filepath.string()));
                throw std::runtime_error(fmt::format("Failed to write downloaded file: {}", filepath.string()));
            }

            progress.completeFile(filename, downloadRes.text.size());
            downloadedFiles.push_back(filepath.string());
            progress.completeOperation();
        }
        // If we have a dataset, download all resources from that dataset
        else if (!M_dataset.empty())
        {
            std::string packageUrl = fmt::format("{}/package_show", M_url);
            
            if (progress.showDebugOutput()) {
                progress.debug("Dataset metadata URL: " + packageUrl);
            }
            
            // Query dataset metadata using CPR
            cpr::Response res = cpr::Get(
                cpr::Url{packageUrl},
                cpr::Parameters{{"id", M_dataset}},
                cpr::Header{{"Authorization", M_apiKey}},
                cpr::VerifySsl{false},
                cpr::Timeout{30000}  // 30 second default timeout
            );

            if (res.status_code != 200)
            {
                progress.error(fmt::format("CKAN dataset query failed: HTTP {}, {}", res.status_code, res.text));
                throw std::runtime_error(fmt::format("CKAN dataset query failed: HTTP {}, {}", res.status_code, res.text));
            }

            // Parse JSON response
            nl::json jsonResponse = nl::json::parse(res.text);
            if (progress.showDebugOutput()) {
                progress.debug("Dataset metadata response: " + jsonResponse.dump(2));
            }
            
            if (!jsonResponse.contains("result") || !jsonResponse["result"].contains("resources"))
            {
                progress.error("Invalid CKAN dataset response structure.");
                throw std::runtime_error("Invalid CKAN dataset response structure.");
            }

            // Download each resource in the dataset
            auto resources = jsonResponse["result"]["resources"];
            int totalFiles = resources.size();
            
            progress.startOperation(fmt::format("Download {} files from CKAN dataset {}", totalFiles, M_dataset));
            
            int fileNum = 0;
            for (const auto& resource : resources)
            {
                fileNum++;
                std::string downloadUrl = resource.value("url", "");
                std::string filename = resource.value("name", "downloaded_resource");
                
                // Handle null or missing size values safely
                std::streamsize fileSize = 0;
                if (resource.contains("size") && !resource["size"].is_null()) {
                    fileSize = resource["size"].get<std::streamsize>();
                }
                
                if (downloadUrl.empty())
                {
                    progress.error("Skipping resource with empty URL: " + filename);
                    continue;
                }

                fs::path filepath = fs::path(dir) / filename;
                
                progress.startFile(filename, fileSize);
                
                if (progress.showDebugOutput()) {
                    progress.debug("Download URL: " + downloadUrl);
                }
                
                // Download the actual file content
                cpr::Response downloadRes = cpr::Get(
                    cpr::Url{downloadUrl},
                    cpr::Header{{"Authorization", M_apiKey}},
                    cpr::VerifySsl{false},
                    cpr::Timeout{30000}  // 30 second default timeout
                );

                if (downloadRes.status_code != 200)
                {
                    progress.error(fmt::format("Failed to download resource {}: HTTP {}", filename, downloadRes.status_code));
                    continue;
                }

                // Write the file content to disk
                std::ofstream file(filepath.string(), std::ios::binary);
                file.write(downloadRes.text.data(), downloadRes.text.size());
                file.close();

                if (!file.good())
                {
                    progress.error(fmt::format("Failed to write file: {}", filepath.string()));
                    continue;
                }

                progress.completeFile(filename, downloadRes.text.size());
                downloadedFiles.push_back(filepath.string());
            }
            
            progress.completeOperation();
        }
        else
        {
            progress.error("No resource ID or dataset specified for CKAN download.");
            throw std::runtime_error("No resource ID or dataset specified for CKAN download.");
        }
    }

    M_worldComm->barrier();
    return downloadedFiles;
}

std::vector<std::string>
RemoteData::CKAN::download(const std::string& dir, int timeout) const
{
    std::vector<std::string> downloadedFiles;
    if (!isInit())
        throw std::runtime_error("CKAN remote data object is not initialized.");

    if (M_worldComm->isMasterRank())
    {
        // Determine progress level
        RemoteDataProgress::Level progressLevel = RemoteDataProgress::Level::NORMAL;
        if (Environment::vm().count("quiet")) {
            progressLevel = RemoteDataProgress::Level::QUIET;
        } else if (Environment::vm().count("debug")) {
            progressLevel = RemoteDataProgress::Level::DEBUG;
        } else if (Environment::vm().count("verbose") || Environment::vm().count("progress")) {
            progressLevel = RemoteDataProgress::Level::VERBOSE;
        }
        
        RemoteDataProgress progress(RemoteDataProgress::Operation::DOWNLOAD, progressLevel);
        
        // If we have a specific resource ID, download that resource
        if (!M_resourceId.empty())
        {
            progress.startOperation("Download CKAN resource " + M_resourceId);
            
            std::string resourceUrl = fmt::format("{}/resource_show", M_url);
            
            auto response = cpr::Get(cpr::Url{resourceUrl},
                                     cpr::Parameters{{"id", M_resourceId}},
                                     cpr::Timeout{timeout});
            
            if (response.status_code != 200)
                throw std::runtime_error(fmt::format("CKAN API request failed: {}", response.status_code));
            
            auto json = nl::json::parse(response.text);
            auto resource = json["result"];
            
            std::string fileName = resource["name"];
            std::string url = resource["url"];
            
            auto progressCallback = [&](cpr::cpr_pf_arg_t downloadTotal, cpr::cpr_pf_arg_t downloadNow, 
                                     cpr::cpr_pf_arg_t uploadTotal, cpr::cpr_pf_arg_t uploadNow, 
                                     intptr_t userdata) -> bool {
                std::streamsize downloadedBytes = static_cast<std::streamsize>(downloadNow);
                std::streamsize totalBytes = static_cast<std::streamsize>(downloadTotal);
                
                static auto lastUpdate = std::chrono::steady_clock::now();
                static cpr::cpr_pf_arg_t lastProgress = 0;
                
                auto now = std::chrono::steady_clock::now();
                auto timeDiff = std::chrono::duration_cast<std::chrono::milliseconds>(now - lastUpdate).count();
                double progressDiff = (totalBytes > 0) ? (static_cast<double>(downloadNow - lastProgress) / totalBytes * 100.0) : 0.0;
                
                bool shouldUpdate = (timeDiff >= 100) || (progressDiff >= 0.5) || (downloadNow == downloadTotal);
                
                if (shouldUpdate) {
                    progress.showProgressBar(fileName, downloadedBytes, totalBytes);
                    lastUpdate = now;
                    lastProgress = downloadNow;
                }
                return true;
            };
            
            auto writeCallback = [&](const std::string_view& data, intptr_t userdata) -> size_t {
                return data.size();
            };
            
            fs::path downloadPath = fs::path(dir) / fileName;
            std::ofstream file(downloadPath, std::ios::binary);
            
            if (!file) {
                throw std::runtime_error(fmt::format("Cannot create file: {}", downloadPath.string()));
            }
            
            auto fileResponse = cpr::Get(cpr::Url{url},
                                         cpr::WriteCallback{writeCallback},
                                         cpr::ProgressCallback{progressCallback},
                                         cpr::Timeout{timeout});
            
            if (fileResponse.status_code == 200) {
                // Redownload to file since we used a custom write callback
                std::ofstream finalFile(downloadPath, std::ios::binary);
                auto finalResponse = cpr::Get(cpr::Url{url}, cpr::Timeout{timeout});
                finalFile << finalResponse.text;
                finalFile.close();
                
                downloadedFiles.push_back(downloadPath.string());
                // progress.info(fmt::format("Downloaded: {} ({} bytes)", fileName, fs::file_size(downloadPath)));
            } else {
                throw std::runtime_error(fmt::format("Download failed: {}", fileResponse.status_code));
            }
            
            progress.completeOperation();
        }
        // Handle dataset downloads (multiple resources)
        else if (!M_dataset.empty())
        {
            progress.startOperation("Download CKAN dataset " + M_dataset);
            
            std::string packageUrl = fmt::format("{}/package_show", M_url);
            
            auto response = cpr::Get(cpr::Url{packageUrl},
                                     cpr::Parameters{{"id", M_dataset}},
                                     cpr::Timeout{timeout});
            
            if (response.status_code != 200)
                throw std::runtime_error(fmt::format("CKAN API request failed: {}", response.status_code));
            
            auto json = nl::json::parse(response.text);
            auto package = json["result"];
            auto resources = package["resources"];
            
            for (const auto& resource : resources)
            {
                std::string fileName = resource["name"];
                std::string url = resource["url"];
                
                auto progressCallback = [&](cpr::cpr_pf_arg_t downloadTotal, cpr::cpr_pf_arg_t downloadNow, 
                                         cpr::cpr_pf_arg_t uploadTotal, cpr::cpr_pf_arg_t uploadNow, 
                                         intptr_t userdata) -> bool {
                    std::streamsize downloadedBytes = static_cast<std::streamsize>(downloadNow);
                    std::streamsize totalBytes = static_cast<std::streamsize>(downloadTotal);
                    
                    static auto lastUpdate = std::chrono::steady_clock::now();
                    static cpr::cpr_pf_arg_t lastProgress = 0;
                    
                    auto now = std::chrono::steady_clock::now();
                    auto timeDiff = std::chrono::duration_cast<std::chrono::milliseconds>(now - lastUpdate).count();
                    double progressDiff = (totalBytes > 0) ? (static_cast<double>(downloadNow - lastProgress) / totalBytes * 100.0) : 0.0;
                    
                    bool shouldUpdate = (timeDiff >= 100) || (progressDiff >= 0.5) || (downloadNow == downloadTotal);
                    
                    if (shouldUpdate) {
                        progress.showProgressBar(fileName, downloadedBytes, totalBytes);
                        lastUpdate = now;
                        lastProgress = downloadNow;
                    }
                    return true;
                };
                
                auto writeCallback = [&](const std::string_view& data, intptr_t userdata) -> size_t {
                    return data.size();
                };
                
                fs::path downloadPath = fs::path(dir) / fileName;
                std::ofstream file(downloadPath, std::ios::binary);
                
                if (!file) {
                    throw std::runtime_error(fmt::format("Cannot create file: {}", downloadPath.string()));
                }
                
                auto fileResponse = cpr::Get(cpr::Url{url},
                                             cpr::WriteCallback{writeCallback},
                                             cpr::ProgressCallback{progressCallback},
                                             cpr::Timeout{timeout});
                
                if (fileResponse.status_code == 200) {
                    // Redownload to file since we used a custom write callback
                    std::ofstream finalFile(downloadPath, std::ios::binary);
                    auto finalResponse = cpr::Get(cpr::Url{url}, cpr::Timeout{timeout});
                    finalFile << finalResponse.text;
                    finalFile.close();
                    
                    downloadedFiles.push_back(downloadPath.string());
                    // progress.info(fmt::format("Downloaded: {} ({} bytes)", fileName, fs::file_size(downloadPath)));
                } else {
                    throw std::runtime_error(fmt::format("Download failed: {}", fileResponse.status_code));
                }
            }
            
            progress.completeOperation();
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

    // Create progress reporter based on command-line options
    RemoteDataProgress::Level progressLevel = RemoteDataProgress::Level::NORMAL;
    if (Environment::vm().count("quiet"))
        progressLevel = RemoteDataProgress::Level::QUIET;
    else if (Environment::vm().count("debug"))
        progressLevel = RemoteDataProgress::Level::DEBUG;
    else if (Environment::vm().count("verbose") || Environment::vm().count("progress"))
        progressLevel = RemoteDataProgress::Level::VERBOSE;
    
    RemoteDataProgress progress(RemoteDataProgress::Operation::UPLOAD, progressLevel);
    progress.startOperation("CKAN", fmt::format("dataset: {}", datasetId));

    // Use default timeout for non-timeout enabled upload method  
    int defaultTimeout = 30000; // 30 seconds
    std::vector<std::string> uploadedResources;
    if (M_worldComm->isMasterRank())
    {
        fs::path dataFsPath(dataPath);
        if (fs::is_regular_file(dataFsPath))
        {
            // Upload single file
            uploadFileWithProgress(dataPath, datasetId, uploadedResources, progress, 1, 1, defaultTimeout);
        }
        else if (fs::is_directory(dataFsPath))
        {
            // Count files for progress reporting
            int fileCount = 0;
            std::vector<fs::directory_entry> files;
            for (const auto& entry : fs::recursive_directory_iterator(dataPath))
            {
                if (entry.is_regular_file())
                {
                    files.push_back(entry);
                    fileCount++;
                }
            }
            
            // Upload each file with progress
            int fileNum = 0;
            for (const auto& entry : files)
            {
                fileNum++;
                uploadFileWithProgress(entry.path().string(), datasetId, uploadedResources, progress, fileNum, fileCount, defaultTimeout);
            }
        }
        else
        {
            progress.error(fmt::format("Unsupported file system object: {}", dataPath));
        }
        
        progress.completeOperation();
    }

    if (sync)
        mpi::broadcast(M_worldComm->globalComm(), uploadedResources, M_worldComm->masterRank());

    return uploadedResources;
}

void
RemoteData::CKAN::uploadFileWithProgress(const std::string& filePath, const std::string& datasetId, std::vector<std::string>& uploadedResources, const RemoteDataProgress& progress, int fileNum, int totalFiles, int timeout) const
{
    fs::path filepath(filePath);
    std::string filename = filepath.filename().string();
    std::streamsize fileSize = fs::file_size(filePath);
    
    progress.startFile(filename, fileSize, fileNum, totalFiles);
    
    std::string url = fmt::format("{}/resource_create", M_url);
    cpr::Header headers = {
        {"Authorization", M_apiKey}
    };

    // Prepare file payload
    std::ifstream fileStream(filePath, std::ios::binary | std::ios::ate);
    fileStream.seekg(0, std::ios::beg);
    std::vector<char> fileData(fileSize);
    fileStream.read(fileData.data(), fileSize);

    // HTTP POST request
    cpr::Multipart multipart = {
        {"package_id", datasetId},
        {"name", filename},
        {"upload", cpr::Buffer(fileData.begin(), fileData.end(), std::move(filename))}
    };

    try
    {
        cpr::Response res = cpr::Post(
            cpr::Url{url},
            headers,
            multipart,
            cpr::VerifySsl{false},
            cpr::Timeout{timeout}
        );
        
        if (res.status_code != 200)
        {
            std::string errorMsg = fmt::format("HTTP {}: {}", res.status_code, res.text.empty() ? "No response" : res.text);
            progress.error(fmt::format("Upload failed for {}: {}", filename, errorMsg));
            return;
        }

        if (res.text.empty())
        {
            progress.error(fmt::format("Upload failed for {}: Empty response from server", filename));
            return;
        }

        // Parse the JSON response for the uploaded resource details
        nl::json jsonResponse = nl::json::parse(res.text);
        
        if (!jsonResponse.contains("result") || !jsonResponse["result"].contains("id"))
        {
            progress.error(fmt::format("Upload failed for {}: Invalid response format", filename));
            return;
        }
        
        std::string resourceId = jsonResponse["result"]["id"];
        uploadedResources.push_back(resourceId);
        
        progress.completeFile(filename, resourceId);
    }
    catch (const std::exception& e)
    {
        progress.error(fmt::format("Upload failed for {}: {}", filename, e.what()));
    }
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
        cpr::Timeout{30000},  // 30 second default timeout
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
        cpr::VerifySsl{false},
        cpr::Timeout{30000}  // 30 second default timeout
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
        cpr::Timeout{30000},  // 30 second default timeout
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
            cpr::Header{{"Authorization", M_apiKey}, {"Content-Type", "application/json"}},
            cpr::Timeout{30000}  // 30 second default timeout
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
                        // Handle null or missing size values safely
                        std::streamsize fileSize = 0;
                        if (resource.contains("size") && !resource["size"].is_null()) {
                            fileSize = resource["size"].get<std::streamsize>();
                        }
                        
                        auto fileInfo = std::make_shared<FileInfo>(
                            resource.value("name", resource.value("id", "unknown")),
                            resource.value("id", ""),
                            fileSize
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

std::vector<std::string>
RemoteData::CKAN::upload(const std::string& dataPath, const std::string& parentId, int timeout) const
{
    CHECK(isInit()) << "CKAN remote data object is not initialized.";
    CHECK(fs::exists(dataPath)) << fmt::format("Data path '{}' does not exist.", dataPath);

    // Use the dataset as the parent ID context for CKAN uploads
    std::string datasetId = parentId.empty() ? M_dataset : parentId;
    
    // Create progress reporter based on command-line options
    RemoteDataProgress::Level progressLevel = RemoteDataProgress::Level::NORMAL;
    if (Environment::vm().count("quiet"))
        progressLevel = RemoteDataProgress::Level::QUIET;
    else if (Environment::vm().count("debug"))
        progressLevel = RemoteDataProgress::Level::DEBUG;
    else if (Environment::vm().count("verbose") || Environment::vm().count("progress"))
        progressLevel = RemoteDataProgress::Level::VERBOSE;
    
    RemoteDataProgress progress(RemoteDataProgress::Operation::UPLOAD, progressLevel);
    progress.startOperation("CKAN", fmt::format("dataset: {}", datasetId));

    std::vector<std::string> uploadedResources;
    if (M_worldComm->isMasterRank())
    {
        fs::path dataFsPath(dataPath);
        if (fs::is_regular_file(dataFsPath))
        {
            // Upload single file
            uploadFileWithProgress(dataPath, datasetId, uploadedResources, progress, 1, 1, timeout);
        }
        else if (fs::is_directory(dataFsPath))
        {
            // Count files for progress reporting
            int fileCount = 0;
            std::vector<fs::directory_entry> files;
            for (const auto& entry : fs::recursive_directory_iterator(dataPath))
            {
                if (entry.is_regular_file())
                {
                    files.push_back(entry);
                    fileCount++;
                }
            }
            
            // Upload each file with progress
            int fileNum = 0;
            for (const auto& entry : files)
            {
                fileNum++;
                uploadFileWithProgress(entry.path().string(), datasetId, uploadedResources, progress, fileNum, fileCount, timeout);
            }
        }
        else
        {
            progress.error(fmt::format("Unsupported file system object: {}", dataPath));
        }
        
        progress.completeOperation();
    }

    return uploadedResources;
}

std::vector<std::string>
RemoteData::CKAN::listOrganizations() const
{
    std::vector<std::string> organizations;
    
    if (!isInit())
    {
        if (M_worldComm->isMasterRank())
            std::cout << "CKAN is not initialized" << std::endl;
        return organizations;
    }

    if (M_worldComm->isMasterRank())
    {
        try 
        {
            // CKAN API endpoint for listing organizations
            std::string url = fmt::format("{}/organization_list", M_url);
            
            // Create headers - API key is optional for public organization list
            cpr::Header headers = {{"Content-Type", "application/json"}};
            if (!M_apiKey.empty()) {
                headers["Authorization"] = M_apiKey;
            }
            
            cpr::Response res = cpr::Get(
                cpr::Url{url},
                headers,
                cpr::VerifySsl{false},
                cpr::Timeout{30000}  // 30 second default timeout
            );

            if (res.status_code == 200)
            {
                nl::json jsonResponse = nl::json::parse(res.text);
                
                if (jsonResponse.contains("success") && jsonResponse["success"].get<bool>())
                {
                    if (jsonResponse.contains("result") && jsonResponse["result"].is_array())
                    {
                        for (const auto& org : jsonResponse["result"])
                        {
                            if (org.is_string())
                            {
                                organizations.push_back(org.get<std::string>());
                            }
                        }
                    }
                }
                else
                {
                    std::string errorMsg = "Unknown error";
                    if (jsonResponse.contains("error"))
                    {
                        errorMsg = jsonResponse["error"].dump();
                    }
                    std::cout << "CKAN organization list failed: " << errorMsg << std::endl;
                }
            }
            else
            {
                std::cout << fmt::format("CKAN organization list failed: HTTP {}, {}", res.status_code, res.text) << std::endl;
            }
        }
        catch (const std::exception& e)
        {
            std::cout << "Error listing CKAN organizations: " << e.what() << std::endl;
        }
    }

    // Broadcast results to all processes
    mpi::broadcast(M_worldComm->globalComm(), organizations, M_worldComm->masterRank());
    
    return organizations;
}

} // namespace Feel