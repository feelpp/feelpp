// Define the test module name
#define BOOST_TEST_MODULE RemoteDataTest


#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/remotedata.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelcore/testsuite.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
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
BOOST_AUTO_TEST_CASE(test_remotedata_girder_delete_if_exist_and_upload)
{
    BOOST_TEST_MESSAGE("Test Girder RemoteData delete if exist and upload");
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

BOOST_AUTO_TEST_SUITE_END()