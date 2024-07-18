#include <feel/feelcore/dataset.hpp>
#include <feel/feelcore/zip.hpp>

#include <feel/feelcore/environment.hpp>

using namespace Feel;

void runSimulation( const Dataset& dataset, const std::string& configDir )
{
    std::cout << fmt::format( "Running dataset: {}", dataset.name ) << std::endl;
    std::cout << fmt::format( "Description: {}", dataset.description ) << std::endl;
    std::cout << fmt::format( "Default configuration: {}", dataset.default_configuration ) << std::endl;
    
    for ( const auto& config : dataset.configurations )
    {
        std::string configPath = configDir + "/" + config.filename;
        if ( fs::exists( configPath ) )
        {
            std::cout << fmt::format( "Running simulation {} with config file: {}", config.name, configPath ) << std::endl;
            //Feel::Environment::setConfigFile( configPath );
            //// Call Feel++ simulation functions here
            //BOOST_TEST_MESSAGE( fmt::format( "Running simulation with config file: {}", configPath ) );
        }
        else
        {
            //BOOST_TEST_MESSAGE( fmt::format( "Config file {} does not exist", configPath ) );
            std::cerr << fmt::format( "Config {} file {} does not exist", config.name, configPath ) << std::endl;
        }
    }
}

int main() 
{
    using namespace Feel;

    std::string zipPath = "/tmp/dataset-fin.zip";
    std::string destDir = "/tmp/dataset-fin";

    extractZipFile(zipPath, destDir);

    std::string jsonPath = destDir + "/fin.json";
    Dataset dataset = loadDatasetDescription(jsonPath);

    runSimulation(dataset, destDir + "/fin");

    return 0;
}