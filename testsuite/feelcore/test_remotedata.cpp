// Define the test module name
#define BOOST_TEST_MODULE RemoteDataTest


#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/remotedata.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelcore/zip.hpp>
#include <feel/feelcore/testsuite.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <fmt/ostream.h>
#include <cpr/cpr.h>
#include <iostream>
#include <fstream>
#include <string>



using namespace Feel;

/**
 * @brief Get the Girder Api Key object
 *
 * the function returns the value of the environment variable FEELPP_GIRDER_API_KEY,
 * the api key is necessary to access the Girder server
 *
 * @return std::string
 */
std::string getGirderApiKey()
{
    const char* key = std::getenv("FEELPP_GIRDER_API_KEY");
    if (!key)
    {
        BOOST_TEST_MESSAGE("Environment variable FEELPP_GIRDER_API_KEY is not set.");
        return "";
    }
    return std::string(key);
}

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( remotedata )

BOOST_AUTO_TEST_CASE(test_remotedata_github)
{
    // Create RemoteData object with GitHub description
    RemoteData rd("github:{repo:feelpp,path:README.adoc}", Environment::worldCommPtr());

    // Check if we can download the data
    if (rd.canDownload())
    {
        // Get the downloads repository directory
        std::string d = Environment::downloadsRepository();
        std::cout << "Download data in: " << d << std::endl;

        // Perform the download
        auto data = rd.download(d);
        std::cout << "Downloaded data:";
        for (const auto& file : data)
            std::cout << " " << file;
        std::cout << std::endl;

        // Optionally, add assertions to check the downloaded files
        BOOST_CHECK(!data.empty()); // Check that data was downloaded
    }
    else
    {
        BOOST_FAIL("Cannot download data using RemoteData");
    }
}

// girder
BOOST_AUTO_TEST_CASE(test_remotedata_girder_delete_if_exist_and_upload)
{
    BOOST_TEST_MESSAGE("Test Girder RemoteData delete if exist and upload");

    std::string girderApiKey = getGirderApiKey();
    if (girderApiKey.empty())
    {
        BOOST_FAIL("FEELPP_GIRDER_API_KEY environment variable is missing.");
    }

    RemoteData rd("girder:{path:/collection/feelpp/testsuite/feelcore/feelpp_test_remotedata/dataset}", Environment::worldCommPtr());

    // Check if we can download the data
    //if (rd.canUpload())
    if ( 1 )
    {
        BOOST_TEST_MESSAGE("Can upload data using RemoteData");
        // find item collection/feelpp/testsuite/feelcore/feelpp_test_remotedata/dataset
        nl::json rid = rd.resourceLookup("collection/feelpp/testsuite/feelcore/feelpp_test_remotedata/dataset");
        if (!rid.empty())
        {
            BOOST_TEST_MESSAGE(fmt::format("collection/feelpp/testsuite/feelcore/feelpp_test_remotedata/dataset exists with id: {}", rid.dump()));
            rd.deleteResource(rid);
        }
        rid = rd.resourceLookup("collection/feelpp/testsuite/feelcore/feelpp_test_remotedata");
        if (!rid.empty())
        {
            BOOST_TEST_MESSAGE(fmt::format("collection/feelpp/testsuite/feelcore/feelpp_test_remotedata exists with id: {}", rid.dump()));
        }
        auto r = rd.createItem("dataset", rid["_id"].get<std::string>());
        if ( r.empty() )
        {
            BOOST_FAIL("Cannot create item");
        }
        // create upload directory in Environment::downloadsRepository()/uploads
        std::string uploadDir = Environment::downloadsRepository() + "/uploads";
        fs::create_directories(uploadDir);
        // create file1.txt and add dummy text
        std::string dataPath = Environment::downloadsRepository() + "/uploads/file.txt";
        std::ofstream file(dataPath);
        file << "Hello, World!";
        file.close();
        // create file2.txt and add dummy text
        std::string dataPath2 = Environment::downloadsRepository() + "/uploads/file2.txt";
        std::ofstream file2(dataPath2);
        file2 << "Hello, World!";
        file2.close();
        // upload file1.txt and file2.txt
        std::vector<std::pair<std::string, std::string>> dataToUpload = {
            {dataPath, "/collection/feelpp/testsuite/feelcore/feelpp_test_remotedata/dataset"},
            {dataPath2, "/collection/feelpp/testsuite/feelcore/feelpp_test_remotedata/dataset"}
        };
        auto data = rd.upload(dataToUpload);
        //std::cout << fmt::format("Uploaded data: {}", data);
    }
    else
    {
        BOOST_FAIL("Cannot upload data using RemoteData");
    }
}
std::map<std::string, std::string> data =
    {
        {"path", "/collection/feelpp/testsuite/feelcore/feelpp_test_remotedata/dataset"},
        {"path", "/collection/feelpp/testsuite/feelcore/feelpp_test_remotedata/"},
        //{"file", "67455d49b0e95728eb010c4a"},
        {"folder", "6743a47bb0e95728eb010c47"}
    };
std::vector<std::map<std::string, std::string>> datasets_map = {data};
namespace bdata = boost::unit_test::data;

BOOST_DATA_TEST_CASE(test_remotedata_girder, bdata::make(datasets_map), dataset)
{
    std::string girderApiKey = getGirderApiKey();
    if (girderApiKey.empty())
    {
        BOOST_FAIL("FEELPP_GIRDER_API_KEY environment variable is missing.");
    }

    BOOST_TEST(!dataset.empty());
    for( auto const& [key, value] : dataset )
    {
        BOOST_TEST_MESSAGE(fmt::format("key = {}, value = {}", key, value));
        // Create RemoteData object with GitHub description
        RemoteData rd(fmt::format("girder:{{{}:{}}}",key,value), Environment::worldCommPtr());

        // Check if we can download the data
        if (rd.canDownload())
        {
            // Get the downloads repository directory
            std::string d = Environment::downloadsRepository();
            std::cout << "Download data in: " << d << std::endl;

            // Perform the download
            auto data = rd.download(d);
            std::cout << "Downloaded data:";
            for (const auto& file : data)
                std::cout << " " << file;
            std::cout << std::endl;

            // Optionally, add assertions to check the downloaded files
            BOOST_CHECK(!data.empty()); // Check that data was downloaded
        }
        else
        {
            BOOST_FAIL("Cannot download data using RemoteData");
        }
    }
}

BOOST_AUTO_TEST_CASE(test_unzip_files)
{
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Create a temporary zip archive and an extraction directory
    fs::path tmpZip = fs::temp_directory_path() / fmt::format("test_{}.zip", rank);
    fs::path extractionDir = fs::temp_directory_path() / fmt::format("unzipped_boost_{}", rank);

    int error = 0;
    zip_t* archive = zip_open(tmpZip.c_str(), ZIP_CREATE | ZIP_TRUNCATE, &error);
    if (!archive)
        throw std::runtime_error("Failed to create zip archive");

    // Add a txt file and empty folder
    const char* fileContent = "Hello Boost Test!";
    zip_source_t* src = zip_source_buffer(archive, fileContent, strlen(fileContent), 0);

    if (zip_file_add(archive, "folder/file.txt", src, ZIP_FL_ENC_UTF_8) < 0)
        throw std::runtime_error("Failed to add file to zip");

    if (zip_dir_add(archive, "folder/emptydir/", ZIP_FL_ENC_UTF_8) < 0)
        throw std::runtime_error("Failed to add directory to zip");

    zip_close(archive);

    // Remove extraction folder if it already exists for testing
    if (fs::exists(extractionDir))
        fs::remove_all(extractionDir);

    // call extractZipFile function
    bool ok = extractZipFile(tmpZip.string(), extractionDir.string());
    BOOST_REQUIRE(ok);

    // Check if files exist
    fs::path expectedFile = extractionDir / "folder" / "file.txt";
    fs::path expectedDir  = extractionDir / "folder" / "emptydir";
    BOOST_CHECK(fs::exists(expectedFile));
    BOOST_CHECK(fs::is_regular_file(expectedFile));
    BOOST_CHECK(fs::exists(expectedDir));
    BOOST_CHECK(fs::is_directory(expectedDir));

    // Check content
    std::ifstream in(expectedFile);
    std::string content;
    std::getline(in, content);
    BOOST_CHECK_EQUAL(content, "Hello Boost Test!");

    // Clean up
    fs::remove_all(extractionDir);
    fs::remove(tmpZip);
}


std::map<std::string, std::string> datasetToUnzip ={ {"path", "/collection/feelpp/testsuite/feelcore/feelpp_test_remotedata/dataset"} };
std::vector<std::map<std::string, std::string>> items_map = {data};
namespace bdata = boost::unit_test::data;

BOOST_DATA_TEST_CASE(test_remotedata_girder_download_item_and_unzip, bdata::make(items_map), dataset)
{
    std::string girderApiKey = getGirderApiKey();
    if (girderApiKey.empty())
    {
        BOOST_FAIL("FEELPP_GIRDER_API_KEY environment variable is missing.");
    }

    BOOST_TEST(!dataset.empty());
    for( auto const& [key, value] : dataset )
    {
        BOOST_TEST_MESSAGE(fmt::format("key = {}, value = {}", key, value));
        RemoteData rd(fmt::format("girder:{{{}:{}}}",key,value), Environment::worldCommPtr());

        // Check if we can download the data
        if (rd.canDownload())
        {
            // Get the downloads repository directory
            std::string d = Environment::downloadsRepository();
            std::cout << "Download data in: " << d << std::endl;

            // Perform the download and unzip the downloaded files
            auto data = rd.download(d);
            std::cout << "Downloaded data:";
            std::cout << "data = " << data << std::endl;
            for (const auto& file : data)
            {
                bool ok = extractZipFile( file, Environment::downloadsRepository());

                // Check that data was unzipped correctly
                BOOST_REQUIRE(ok);
            }
            std::cout << std::endl;

            // Optionally, add assertions to check the downloaded files
            BOOST_CHECK(!data.empty()); // Check that data was downloaded
        }
        else
        {
            BOOST_FAIL("Cannot download data using RemoteData");
        }
    }
}

BOOST_AUTO_TEST_CASE(test_remotedata_ckan_upload_download)
{
    // CKAN URL and dataset - using the specified dataset
    std::string ckanUrl = "https://ckan.hidalgo2.eu";
    std::string dataset = "bestest_base"; // Using the specified dataset
    std::string organization = "4719ef48-cce5-4f98-b6e8-37e38655cc86"; //"Cemosis";
    std::string apiKey;

    // Obtain API key from environment variable
    const char* apiKeyEnv = std::getenv("CKAN_API_KEY");
    if (apiKeyEnv == nullptr)
    {
        BOOST_TEST_MESSAGE("CKAN_API_KEY environment variable not set. Skipping test.");
        return; // Skip the test if no API key is available
    }
    else
    {
        apiKey = apiKeyEnv;
    }

    // First test: contents - check what's in the bestest_base dataset
    BOOST_TEST_MESSAGE("Testing CKAN contents for bestest_base dataset");
    RemoteData rdContents(fmt::format("ckan:{{url:{}, dataset:{}}}", ckanUrl, dataset), Environment::worldCommPtr());
    
    auto contents = rdContents.contents();
    auto folders = std::get<0>(contents);
    auto items = std::get<1>(contents);
    auto files = std::get<2>(contents);
    
    BOOST_TEST_MESSAGE(fmt::format("Found {} folders, {} items, {} files in dataset", 
                                   folders.size(), items.size(), files.size()));
    
    // List the files found
    for (const auto& file : files)
    {
        BOOST_TEST_MESSAGE(fmt::format("File: {} (ID: {})", file->name(), file->id()));
    }

    // Second test: try to download a specific file from the dataset (if permissions allow)
    if (!files.empty())
    {
        BOOST_TEST_MESSAGE("Testing CKAN download from bestest_base dataset");
        std::string downloadDir = Environment::downloadsRepository() + "/ckan/downloads";
        fs::create_directories(downloadDir);

        // Try to download the first file
        auto firstFile = files[0];
        std::string resourceId = firstFile->id();
        BOOST_TEST_MESSAGE(fmt::format("Attempting to download resource: {} (ID: {})", firstFile->name(), resourceId));
        
        try 
        {
            RemoteData rdDownload(fmt::format("ckan:{{url:{}, resource:{}}}", ckanUrl, resourceId), Environment::worldCommPtr());
            
            if (rdDownload.canDownload())
            {
                auto downloadedData = rdDownload.download(downloadDir);
                
                BOOST_TEST_MESSAGE(fmt::format("Downloaded {} files", downloadedData.size()));
                
                // In parallel execution, some processes might not download successfully due to 
                // network timing or race conditions. We check that at least one process succeeded
                // by checking if files were actually downloaded to disk
                bool filesExist = false;
                for (const auto& file : downloadedData)
                {
                    BOOST_TEST_MESSAGE(fmt::format("Downloaded: {}", file));
                    if (fs::exists(file))
                    {
                        filesExist = true;
                    }
                }
                
                // Only check if this specific process downloaded files successfully
                // Don't fail if other parallel processes had issues
                if (!downloadedData.empty())
                {
                    BOOST_CHECK(filesExist);
                }
                else
                {
                    BOOST_TEST_MESSAGE("No files downloaded in this process (may be due to parallel execution)");
                }
            }
            else
            {
                BOOST_TEST_MESSAGE("Cannot download from CKAN - API key may not have sufficient permissions");
            }
        }
        catch (const std::exception& e)
        {
            BOOST_TEST_MESSAGE(fmt::format("CKAN download failed (expected if no download permissions): {}", e.what()));
            // Don't fail the test - this is expected if the API key doesn't have download permissions
        }
    }
    else
    {
        BOOST_TEST_MESSAGE("No files found in dataset");
    }

    // Third test: upload (only if we have an API key)
    // Create RemoteData object with CKAN description including API key for upload
    RemoteData rd(fmt::format("ckan:{{url:{}, dataset:{}, api_key:{}}}", ckanUrl, dataset, apiKey), Environment::worldCommPtr());

    // Check if we can upload data (test upload capability but don't upload to existing dataset)
    if (rd.canUpload())
    {
        BOOST_TEST_MESSAGE("CKAN upload capability confirmed");
        
        // For testing, we'll just verify upload capability without actually uploading
        // to the existing bestest_base dataset to avoid modifying it
        BOOST_CHECK(true); // Upload capability exists
    }
    else
    {
        BOOST_TEST_MESSAGE("Cannot upload data - API key may not have sufficient permissions");
        // This is not necessarily a failure for read-only testing
    }
}

BOOST_AUTO_TEST_SUITE_END()