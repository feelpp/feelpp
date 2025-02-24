#include <feel/feelcore/json.hpp>
#include <fstream>

namespace Feel
{
using json = nlohmann::json;

struct Configuration {
    std::string name;
    std::string filename;
    std::string description;
};

struct Dataset {
    std::string name;
    std::string description;
    std::string default_configuration;
    std::vector<Configuration> configurations;
    std::map<std::string, std::string> metadata;
};

Dataset loadDatasetDescription(const std::string& jsonPath) 
{
    std::ifstream ifs(jsonPath);
    json j;
    ifs >> j;

    Dataset dataset;
    dataset.name = j.at("name").get<std::string>();
    dataset.default_configuration = j.at( "default" ).get<std::string>();
    dataset.description = j.at("description").get<std::string>();
    for (const auto& config : j.at("configurations")) {
        dataset.configurations.push_back({
            config.at("name").get<std::string>(),
            config.at("filename").get<std::string>(),
            config.at("description").get<std::string>()
        });
    }
    dataset.metadata = j.at("metadata").get<std::map<std::string, std::string>>();
    return dataset;
}
} // namespace Feel