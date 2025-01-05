// Define the test module name
#define BOOST_TEST_MODULE RemoteDataTest

#include <feel/feelcore/testsuite.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/remotedata.hpp>
#include <feel/feelcore/json.hpp>
#include <cpr/cpr.h>
#include <iostream>
#include <fstream>
#include <string>



using namespace Feel;

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
BOOST_AUTO_TEST_CASE(test_remotedata_girder)
{

    // Girder API base URL
    std::string girderApiUrl = "https://girder.math.unistra.fr/api/v1";

    // The path to the file in Girder
    std::string filePath = "/collection/feelpp/testsuite/feelcore/feelpp_test_remotedata/dataset";

    // If authentication is required, set your token
    // std::string authToken = "your_token_here";

    // 1. Lookup the resource to get the file ID
    cpr::Response res = cpr::Get(
        cpr::Url{girderApiUrl + "/resource/lookup"},
        cpr::Parameters{{"path", filePath}}
        // If authentication is needed:
        // , cpr::Header{{"Girder-Token", authToken}}
    );

    if (res.status_code != 200) 
    {
        BOOST_FAIL( fmt::format("Failed to lookup resource. HTTP status code: {}\n error message: {} ", res.status_code, res.text) );
    }

    // Parse the JSON response to get the file ID
    nlohmann::json jsonResponse = nlohmann::json::parse(res.text);
    if (jsonResponse.contains("_id")) 
    {
        std::string fileId = jsonResponse["_id"];
        std::string modelType = jsonResponse["_modelType"];
        std::string name = jsonResponse["name"];
        std::cout << fmt::format("File ID: {}, Model Type: {} ", fileId, modelType) << std::endl;
        std::string url =  fmt::format("{}/{}/{}/download",girderApiUrl,modelType,fileId);
        std::cout << "Download URL: " << url << std::endl;
        // 2. Download the file using the file ID
        cpr::Response fileRes = cpr::Get(
            cpr::Url{url},
            // If authentication is needed:
            // cpr::Header{{"Girder-Token", authToken}}
            cpr::VerifySsl{false} // Add this if SSL verification causes issues
        );

        if (fileRes.status_code != 200) 
        {
            BOOST_FAIL( fmt::format("Failed to download file. HTTP status code: {}\n error message: {} ", fileRes.status_code, fileRes.text) );
        }

        // Save the file content to a local file
        std::string outputFileName = fmt::format("{}.zip",name); // Set desired output file name
        std::ofstream outputFile(outputFileName, std::ios::binary);
        outputFile << fileRes.text;
        outputFile.close();

        std::cout << "File downloaded successfully: " << outputFileName << std::endl;
    } else {
        std::cerr << "File ID not found in the response." << std::endl;
        return ;
    }


}

BOOST_AUTO_TEST_SUITE_END()