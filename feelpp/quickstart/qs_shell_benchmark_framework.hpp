//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
#ifndef FEELPP_QUICKSTART_QS_SHELL_BENCHMARK_FRAMEWORK_HPP
#define FEELPP_QUICKSTART_QS_SHELL_BENCHMARK_FRAMEWORK_HPP 1

#include <feel/feelcore/json.hpp>
#include <feel/feeldiscr/mesh_fwd.hpp>
#include <feel/feelmesh/hypercube.hpp>

#include <Eigen/Core>

#include <array>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace Feel::Quickstart::ShellBenchmark
{
using json = nl::json;
using mesh_type = Mesh<Hypercube<3>>;
using mesh_ptrtype = std::shared_ptr<mesh_type>;
using shell_vec = Eigen::Matrix<double, 3, 1>;

enum class GeometryKind
{
    Rectangular,
    Circular
};

struct BenchmarkConfig
{
    std::string name;
    std::string description;
    GeometryKind geometryKind = GeometryKind::Rectangular;
    double length = 1.0;
    double width = 1.0;
    double radius = 0.5;
    double thickness = 0.1;
    std::string meshTemplatePath;
    int nx = 1;
    int ny = 1;
    int nz = 1;
    int nr = 1;
    int nt = 1;
    double E = 1.0e6;
    double nu = 0.3;
    shell_vec probe = shell_vec::Zero();
    std::string probeLabel = "probe";
    shell_vec probeDirection = shell_vec::UnitZ();
    std::string bodyForceExpression = "{0,0,0}";
    std::vector<std::string> clampMarkers;

    struct PointConstraint
    {
        std::string marker;
        std::array<bool, 3> components = { true, true, true };
        shell_vec value = shell_vec::Zero();
    };

    struct FacePressureLoad
    {
        std::string marker;
        double value = 0.0;
    };

    struct FaceTractionLoad
    {
        std::string marker;
        std::string expression = "{0,0,0}";
    };

    struct FaceTotalForceLoad
    {
        std::string marker;
        shell_vec value = shell_vec::Zero();
    };

    struct PointLoad
    {
        std::string marker;
        shell_vec value = shell_vec::Zero();
    };

    struct ReferenceValue
    {
        bool hasProbeValue = false;
        double probeValue = 0.0;
        bool hasAbsoluteTolerance = false;
        double absoluteTolerance = 0.0;
        bool hasRelativeTolerance = false;
        double relativeTolerance = 0.0;
        std::string description;
    };

    std::vector<FacePressureLoad> pressureLoads;
    std::vector<FaceTractionLoad> tractionLoads;
    std::vector<FaceTotalForceLoad> totalForceLoads;
    std::vector<PointConstraint> pointConstraints;
    std::vector<PointLoad> pointLoads;
    ReferenceValue reference;
};

struct MeshBuildResult
{
    mesh_ptrtype mesh;
    shell_vec probe;
};

struct MeshBuildOptions
{
    std::string repositoryVariantPath;
    std::string meshVariantTag;
    bool requireSingleLayer = true;
};

template<typename T>
T
jsonValue( json const& object, std::string const& key, T const& defaultValue )
{
    if ( !object.contains( key ) )
        return defaultValue;

    return object.at( key ).get<T>();
}

json const& jsonArray( json const& object, std::string const& key );
std::string shellVectorExpression( shell_vec const& value );
json loadBenchmarkSpecs( std::string const& path );
json const& benchmarkSpec( json const& specs, std::string const& name );
BenchmarkConfig benchmarkPreset( json const& specs, std::string const& name );
MeshBuildResult buildMesh( BenchmarkConfig const& cfg, MeshBuildOptions const& options = {} );
} // namespace Feel::Quickstart::ShellBenchmark

#endif
