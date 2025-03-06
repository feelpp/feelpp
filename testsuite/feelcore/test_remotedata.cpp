// Define the test module name
#define BOOST_TEST_MODULE RemoteDataTest


#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/remotedata.hpp>
#include <feel/feelcore/json.hpp>
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

BOOST_AUTO_TEST_CASE(test_remotedata_ckan_upload_download)
{
    // CKAN URL and organization
    std::string ckanUrl = "https://ckan.hidalgo2.eu";
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

    // Create RemoteData object with CKAN description
    RemoteData rd(fmt::format("ckan:{{url:{}, organization:{}, api_key:{}}}", ckanUrl, organization, apiKey), Environment::worldCommPtr());

    // Check if we can upload data
    if (rd.canUpload())
    {
        // Prepare data to upload
        std::string uploadDir = Environment::downloadsRepository() + "/ckan/uploads";
        fs::create_directories(uploadDir);

        // Create a test file to upload
        std::string dataPath = uploadDir + "/file.txt";
        std::ofstream file(dataPath);
        file << "Hello, CKAN!";
        file.close();

        // Create a unique dataset name
        std::string datasetName = fmt::format("test_dataset_{}", Environment::worldComm().rank());
        nl::json dataset = rd.createDataset(datasetName);
        if (dataset.empty())
        {
            BOOST_FAIL("Cannot create dataset on CKAN");
        }

        std::string datasetId = dataset["id"].get<std::string>();
        BOOST_TEST_MESSAGE(fmt::format("Created dataset with ID: {}", datasetId));
        // Upload the file to the dataset
        auto uploadedResources = rd.upload(dataPath, datasetId, true); // sync = true
        //std::cout << fmt::format("Uploaded resources: {}", uploadedResources) << std::endl;
        BOOST_CHECK(!uploadedResources.empty());

        // Now download the data
        if (rd.canDownload())
        {
            std::string downloadDir = Environment::downloadsRepository() + "/ckan/";
            fs::create_directories(downloadDir);

            // Create a RemoteData object for downloading, specifying the dataset ID
            RemoteData rdDownload(fmt::format("ckan:{{url:{}, organization: {} datasetId:{}}}", ckanUrl, organization, datasetId), Environment::worldCommPtr());

            auto downloadedData = rdDownload.download(downloadDir);
            BOOST_CHECK(!downloadedData.empty());

            // Verify the downloaded file
            std::string downloadedFilePath = downloadDir + "/file.txt";
            BOOST_CHECK(fs::exists(downloadedFilePath));

            // Read the content and verify
            std::ifstream downloadedFile(downloadedFilePath);
            std::string content((std::istreambuf_iterator<char>(downloadedFile)), std::istreambuf_iterator<char>());
            BOOST_CHECK_EQUAL(content, "Hello, CKAN!");
        }
        else
        {
            BOOST_FAIL("Cannot download data using RemoteData");
        }

        // Clean up by deleting the dataset
        bool deleteSuccess = rd.deleteDataset(datasetId);
        BOOST_CHECK(deleteSuccess);
    }
    else
    {
        BOOST_FAIL("Cannot upload data using RemoteData");
    }
}

BOOST_AUTO_TEST_SUITE_END()