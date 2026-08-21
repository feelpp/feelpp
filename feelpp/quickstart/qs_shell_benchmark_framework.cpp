/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#include "qs_shell_benchmark_framework.hpp"

#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/geo.hpp>

#include <boost/format.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numbers>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace Feel::Quickstart::ShellBenchmark
{
namespace
{
json const&
jsonObject( json const& object, std::string const& key )
{
    if ( !object.contains( key ) || !object.at( key ).is_object() )
        throw std::invalid_argument( "missing json object '" + key + "'" );

    return object.at( key );
}

std::string
formatDouble( double value )
{
    std::ostringstream os;
    os << std::setprecision( 16 ) << value;
    return os.str();
}

std::string
jsonExpression( json const& object, std::string const& key, std::string const& defaultValue )
{
    if ( !object.contains( key ) )
        return defaultValue;

    auto const& value = object.at( key );
    if ( value.is_string() )
        return value.get<std::string>();
    if ( value.is_number() )
        return formatDouble( value.get<double>() );

    throw std::invalid_argument( "expected '" + key + "' to be a scalar expression" );
}

shell_vec
jsonShellVector( json const& value, std::string const& description )
{
    if ( !value.is_array() || value.size() != 3 )
        throw std::invalid_argument( description + " must be a size-3 array" );

    shell_vec v;
    for ( int i = 0; i < 3; ++i )
        v( i ) = value.at( i ).get<double>();
    return v;
}

shell_vec
jsonShellVector( json const& object, std::string const& key, shell_vec const& defaultValue )
{
    if ( !object.contains( key ) )
        return defaultValue;

    return jsonShellVector( object.at( key ), key );
}

shell_voigt
jsonShellVoigt( json const& value, std::string const& description )
{
    if ( !value.is_array() || value.size() != 6 )
        throw std::invalid_argument( description + " must be a size-6 array" );

    shell_voigt v;
    for ( int i = 0; i < 6; ++i )
        v( i ) = value.at( i ).get<double>();
    return v;
}

std::array<bool, 3>
jsonComponentMask( json const& object, std::string const& key, std::array<bool, 3> const& defaultValue )
{
    if ( !object.contains( key ) )
        return defaultValue;

    auto const& value = object.at( key );
    if ( !value.is_array() || value.size() != 3 )
        throw std::invalid_argument( key + " must be a size-3 array" );

    std::array<bool, 3> mask = defaultValue;
    for ( int i = 0; i < 3; ++i )
        mask[i] = value.at( i ).get<bool>();
    return mask;
}

GeometryKind
geometryKindFromString( std::string const& kind )
{
    if ( kind == "rectangular" )
        return GeometryKind::Rectangular;
    if ( kind == "circular" )
        return GeometryKind::Circular;

    throw std::invalid_argument( "unsupported geometry kind '" + kind + "'" );
}

std::string
defaultMeshTemplatePath( GeometryKind kind )
{
    switch ( kind )
    {
    case GeometryKind::Rectangular:
        return "$top_srcdir/feelpp/quickstart/meshes/sb9_rectangular_shell.geo.tmpl";
    case GeometryKind::Circular:
        return "$top_srcdir/feelpp/quickstart/meshes/sb9_circular_shell.geo.tmpl";
    }

    throw std::invalid_argument( "unsupported geometry kind for default mesh template" );
}

std::string
loadTextFile( std::string const& path )
{
    std::ifstream input( path );
    if ( !input )
        throw std::invalid_argument( "unable to open geo template '" + path + "'" );

    std::ostringstream content;
    content << input.rdbuf();
    return content.str();
}

std::string
replaceAll( std::string text, std::string const& from, std::string const& to )
{
    std::size_t pos = 0;
    while ( ( pos = text.find( from, pos ) ) != std::string::npos )
    {
        text.replace( pos, from.size(), to );
        pos += to.size();
    }
    return text;
}

std::string
makeGeoDescription( BenchmarkConfig const& cfg )
{
    auto geoDesc = loadTextFile( Environment::expand( cfg.meshTemplatePath ) );
    double const z0 = -0.5 * cfg.thickness;
    double const z1 = 0.5 * cfg.thickness;
    double const innerRadius = cfg.radius / 3.0;

    std::array<std::pair<char const*, std::string>, 11> replacements = {{
        { "@LENGTH@", formatDouble( cfg.length ) },
        { "@WIDTH@", formatDouble( cfg.width ) },
        { "@RADIUS@", formatDouble( cfg.radius ) },
        { "@INNER_RADIUS@", formatDouble( innerRadius ) },
        { "@THICKNESS@", formatDouble( cfg.thickness ) },
        { "@Z0@", formatDouble( z0 ) },
        { "@Z1@", formatDouble( z1 ) },
        { "@NX@", std::to_string( cfg.nx + 1 ) },
        { "@NY@", std::to_string( cfg.ny + 1 ) },
        { "@NZ@", std::to_string( cfg.nz ) },
        { "@NR@", std::to_string( cfg.nr + 1 ) }
    }};
    for ( auto const& [token, value] : replacements )
        geoDesc = replaceAll( geoDesc, token, value );

    geoDesc = replaceAll( geoDesc, "@NT@", std::to_string( cfg.nt + 1 ) );
    return geoDesc;
}

std::filesystem::path
resolveSpecPath( std::string const& includePath, std::string const& parentPath )
{
    std::filesystem::path resolved( Environment::expand( includePath ) );
    if ( resolved.is_absolute() )
        return resolved.lexically_normal();

    std::filesystem::path parent( parentPath );
    if ( std::filesystem::is_regular_file( parent ) )
        parent = parent.parent_path();
    return ( parent / resolved ).lexically_normal();
}

void
mergeSpecEntries( json const& source,
                  std::string const& key,
                  json& destination,
                  std::string const& sourcePath )
{
    if ( !source.contains( key ) )
        return;

    auto const& entries = jsonObject( source, key );
    for ( auto const& [name, spec] : entries.items() )
    {
        if ( destination.contains( name ) )
            throw std::invalid_argument( "duplicate benchmark/validation '" + name + "' while merging '" + sourcePath + "'" );
        destination[name] = spec;
    }
}

json
loadBenchmarkSpecsRecursive( std::string const& path, std::set<std::string>& visited )
{
    if ( !visited.insert( path ).second )
        throw std::invalid_argument( "cyclic benchmark specification include detected at '" + path + "'" );

    std::ifstream specsStream( path );
    if ( !specsStream )
        throw std::invalid_argument( "unable to open benchmark specs '" + path + "'" );

    auto root = json::parse( specsStream );
    json merged = json::object();
    merged["benchmarks"] = json::object();

    if ( root.contains( "include" ) )
    {
        for ( auto const& includeEntry : jsonArray( root, "include" ) )
        {
            if ( !includeEntry.is_string() )
                throw std::invalid_argument( "expected include entries to be strings in '" + path + "'" );

            auto const includePath = resolveSpecPath( includeEntry.get<std::string>(), path );
            auto const child = loadBenchmarkSpecsRecursive( includePath.string(), visited );
            mergeSpecEntries( child, "benchmarks", merged["benchmarks"], includePath.string() );
        }
    }

    mergeSpecEntries( root, "benchmarks", merged["benchmarks"], path );
    mergeSpecEntries( root, "validations", merged["benchmarks"], path );
    return merged;
}

std::optional<double>
benchmarkHsizeOverride()
{
    constexpr double feelppDefaultBenchmarkH = 0.1;
    auto const& vm = Environment::vm();
    auto const it = vm.find( "benchmark.hsize" );
    if ( it == vm.end() || it->second.defaulted() )
        return std::nullopt;

    double const requestedH = it->second.as<double>();
    if ( requestedH <= 0.0 )
        return std::nullopt;
    if ( std::abs( requestedH - feelppDefaultBenchmarkH ) <= 1.0e-12 )
        return std::nullopt;

    return requestedH;
}

std::string
normalizePathFragment( std::string fragment )
{
    std::replace( fragment.begin(), fragment.end(), '\\', '/' );
    while ( !fragment.empty() && fragment.front() == '/' )
        fragment.erase( fragment.begin() );
    while ( !fragment.empty() && fragment.back() == '/' )
        fragment.pop_back();
    return fragment;
}
} // namespace

json const&
jsonArray( json const& object, std::string const& key )
{
    if ( !object.contains( key ) || !object.at( key ).is_array() )
        throw std::invalid_argument( "missing json array '" + key + "'" );

    return object.at( key );
}

std::string
shellVectorExpression( shell_vec const& value )
{
    std::ostringstream os;
    os << "{"
       << formatDouble( value( 0 ) ) << ","
       << formatDouble( value( 1 ) ) << ","
       << formatDouble( value( 2 ) ) << "}";
    return os.str();
}

json
loadBenchmarkSpecs( std::string const& path )
{
    std::set<std::string> visited;
    return loadBenchmarkSpecsRecursive( resolveSpecPath( path, std::filesystem::current_path().string() ).string(), visited );
}

json const&
benchmarkSpec( json const& specs, std::string const& name )
{
    return jsonObject( jsonObject( specs, "benchmarks" ), name );
}

BenchmarkConfig
benchmarkPreset( json const& specs, std::string const& name )
{
    auto const& benchmarkSpecs = benchmarkSpec( specs, name );
    auto const& geometrySpecs = jsonObject( benchmarkSpecs, "geometry" );
    auto const& materialSpecs = jsonObject( benchmarkSpecs, "material" );

    BenchmarkConfig cfg;
    cfg.name = name;
    cfg.description = jsonValue<std::string>( benchmarkSpecs, "description", "" );
    cfg.geometryKind = geometryKindFromString( jsonValue<std::string>( geometrySpecs, "kind", "rectangular" ) );
    cfg.meshTemplatePath = jsonValue<std::string>( geometrySpecs, "template", defaultMeshTemplatePath( cfg.geometryKind ) );
    cfg.length = jsonValue<double>( geometrySpecs, "length", cfg.length );
    cfg.width = jsonValue<double>( geometrySpecs, "width", cfg.width );
    cfg.radius = jsonValue<double>( geometrySpecs, "radius", cfg.radius );
    cfg.thickness = jsonValue<double>( geometrySpecs, "thickness", cfg.thickness );
    cfg.nx = jsonValue<int>( geometrySpecs, "nx", cfg.nx );
    cfg.ny = jsonValue<int>( geometrySpecs, "ny", cfg.ny );
    cfg.nz = jsonValue<int>( geometrySpecs, "nz", cfg.nz );
    cfg.nr = jsonValue<int>( geometrySpecs, "nr", cfg.nr );
    cfg.nt = jsonValue<int>( geometrySpecs, "nt", cfg.nt );
    cfg.E = jsonValue<double>( materialSpecs, "E", cfg.E );
    cfg.nu = jsonValue<double>( materialSpecs, "nu", cfg.nu );

    if ( benchmarkSpecs.contains( "probe" ) )
    {
        auto const& probeSpecs = jsonObject( benchmarkSpecs, "probe" );
        cfg.probe = jsonShellVector( probeSpecs, "point", cfg.probe );
        cfg.probeLabel = jsonValue<std::string>( probeSpecs, "label", cfg.probeLabel );
        cfg.probeDirection = jsonShellVector( probeSpecs, "direction", cfg.probeDirection );
    }

    if ( benchmarkSpecs.contains( "boundaryConditions" ) )
    {
        auto const& bcSpecs = jsonObject( benchmarkSpecs, "boundaryConditions" );
        if ( bcSpecs.contains( "clampMarkers" ) )
            cfg.clampMarkers = bcSpecs.at( "clampMarkers" ).get<std::vector<std::string>>();

        if ( bcSpecs.contains( "pointConstraints" ) )
        {
            for ( auto const& constraintSpec : bcSpecs.at( "pointConstraints" ) )
            {
                cfg.pointConstraints.push_back(
                    { constraintSpec.at( "marker" ).get<std::string>(),
                      jsonComponentMask( constraintSpec, "components", { true, true, true } ),
                      jsonShellVector( constraintSpec, "value", shell_vec::Zero() ) } );
            }
        }
    }

    if ( benchmarkSpecs.contains( "loads" ) )
    {
        auto const& loadSpecs = jsonObject( benchmarkSpecs, "loads" );
        cfg.bodyForceExpression = jsonExpression( loadSpecs, "bodyForce", cfg.bodyForceExpression );

        if ( loadSpecs.contains( "pressures" ) )
        {
            for ( auto const& pressureSpec : loadSpecs.at( "pressures" ) )
            {
                cfg.pressureLoads.push_back( { pressureSpec.at( "marker" ).get<std::string>(),
                                               pressureSpec.at( "value" ).get<double>() } );
            }
        }

        if ( loadSpecs.contains( "tractions" ) )
        {
            for ( auto const& tractionSpec : loadSpecs.at( "tractions" ) )
            {
                cfg.tractionLoads.push_back( { tractionSpec.at( "marker" ).get<std::string>(),
                                               jsonExpression( tractionSpec, "expr", "{0,0,0}" ) } );
            }
        }

        if ( loadSpecs.contains( "totalForces" ) )
        {
            for ( auto const& totalForceSpec : loadSpecs.at( "totalForces" ) )
            {
                cfg.totalForceLoads.push_back( { totalForceSpec.at( "marker" ).get<std::string>(),
                                                 jsonShellVector( totalForceSpec.at( "value" ), "total force value" ) } );
            }
        }

        if ( loadSpecs.contains( "pointLoads" ) )
        {
            for ( auto const& pointLoadSpec : loadSpecs.at( "pointLoads" ) )
            {
                cfg.pointLoads.push_back( { pointLoadSpec.at( "marker" ).get<std::string>(),
                                            jsonShellVector( pointLoadSpec.at( "value" ), "point load value" ),
                                            jsonValue<std::string>( pointLoadSpec, "quantity", "per-point" ) } );
            }
        }

        if ( loadSpecs.contains( "pointMoments" ) )
        {
            for ( auto const& pointMomentSpec : loadSpecs.at( "pointMoments" ) )
            {
                cfg.pointMoments.push_back( { pointMomentSpec.at( "marker" ).get<std::string>(),
                                              jsonShellVector( pointMomentSpec.at( "value" ), "point moment value" ),
                                              jsonValue<std::string>( pointMomentSpec, "quantity", "per-point" ) } );
            }
        }
    }

    if ( benchmarkSpecs.contains( "reference" ) )
    {
        auto const& referenceSpecs = jsonObject( benchmarkSpecs, "reference" );
        if ( referenceSpecs.contains( "probeValue" ) )
        {
            cfg.reference.hasProbeValue = true;
            cfg.reference.probeValue = referenceSpecs.at( "probeValue" ).get<double>();
        }
        if ( referenceSpecs.contains( "absoluteTolerance" ) )
        {
            cfg.reference.hasAbsoluteTolerance = true;
            cfg.reference.absoluteTolerance = referenceSpecs.at( "absoluteTolerance" ).get<double>();
        }
        if ( referenceSpecs.contains( "relativeTolerance" ) )
        {
            cfg.reference.hasRelativeTolerance = true;
            cfg.reference.relativeTolerance = referenceSpecs.at( "relativeTolerance" ).get<double>();
        }
        cfg.reference.description = jsonValue<std::string>( referenceSpecs, "description", "" );
    }

    if ( benchmarkSpecs.contains( "fieldReferences" ) )
    {
        for ( auto const& fieldReferenceSpec : benchmarkSpecs.at( "fieldReferences" ) )
        {
            BenchmarkConfig::FieldReference fieldReference;
            fieldReference.name = jsonValue<std::string>( fieldReferenceSpec, "name", "" );
            fieldReference.description = jsonValue<std::string>( fieldReferenceSpec, "description", "" );
            fieldReference.point = jsonShellVector( fieldReferenceSpec, "point", fieldReference.point );
            if ( fieldReferenceSpec.contains( "absoluteTolerance" ) )
            {
                fieldReference.hasAbsoluteTolerance = true;
                fieldReference.absoluteTolerance = fieldReferenceSpec.at( "absoluteTolerance" ).get<double>();
            }
            if ( fieldReferenceSpec.contains( "epsilonAbsoluteTolerances" ) )
            {
                fieldReference.hasEpsilonAbsoluteTolerances = true;
                fieldReference.epsilonAbsoluteTolerances =
                    jsonShellVoigt( fieldReferenceSpec.at( "epsilonAbsoluteTolerances" ),
                                    "field reference epsilon absolute tolerances" );
            }
            if ( fieldReferenceSpec.contains( "sigmaAbsoluteTolerances" ) )
            {
                fieldReference.hasSigmaAbsoluteTolerances = true;
                fieldReference.sigmaAbsoluteTolerances =
                    jsonShellVoigt( fieldReferenceSpec.at( "sigmaAbsoluteTolerances" ),
                                    "field reference sigma absolute tolerances" );
            }
            if ( fieldReferenceSpec.contains( "relativeTolerance" ) )
            {
                fieldReference.hasRelativeTolerance = true;
                fieldReference.relativeTolerance = fieldReferenceSpec.at( "relativeTolerance" ).get<double>();
            }

            auto const& targetsSpec = jsonObject( fieldReferenceSpec, "targets" );
            for ( auto const& [targetName, targetSpec] : targetsSpec.items() )
            {
                fieldReference.targets.emplace(
                    targetName,
                    BenchmarkConfig::FieldReferenceTarget{
                        jsonShellVoigt( targetSpec.at( "epsilon" ), "field reference epsilon" ),
                        jsonShellVoigt( targetSpec.at( "sigma" ), "field reference sigma" ) } );
            }
            cfg.fieldReferences.push_back( std::move( fieldReference ) );
        }
    }

    if ( cfg.clampMarkers.empty() && cfg.pointConstraints.empty() )
        throw std::invalid_argument( "benchmark '" + name + "' has no displacement constraints" );

    return cfg;
}

MeshBuildResult
buildMesh( BenchmarkConfig const& cfg, MeshBuildOptions const& options )
{
    auto meshCfg = cfg;
    if ( options.requireSingleLayer && meshCfg.nz != 1 )
    {
        throw std::invalid_argument( Environment::about().appName() +
                                     " expects exactly one hexahedral layer through the thickness" );
    }

    auto const requestedH = benchmarkHsizeOverride();
    if ( requestedH )
    {
        if ( meshCfg.geometryKind == GeometryKind::Rectangular )
        {
            meshCfg.nx = std::max( 1, static_cast<int>( std::ceil( meshCfg.length / *requestedH ) ) );
            meshCfg.ny = std::max( 1, static_cast<int>( std::ceil( meshCfg.width / *requestedH ) ) );
        }
        else
        {
            double const innerRadius = meshCfg.radius / 3.0;
            double const radialSpan = std::max( meshCfg.radius - innerRadius, *requestedH );
            double const quarterArcLength = 0.5 * std::numbers::pi_v<double> * meshCfg.radius;
            meshCfg.nr = std::max( 1, static_cast<int>( std::ceil( radialSpan / *requestedH ) ) );
            meshCfg.nt = std::max( 1, static_cast<int>( std::ceil( quarterArcLength / *requestedH ) ) );
        }
    }

    std::ostringstream repository;
    repository << "quickstart/" << Environment::about().appName() << "/" << meshCfg.name;
    auto const repositoryVariantPath = normalizePathFragment( options.repositoryVariantPath );
    if ( !repositoryVariantPath.empty() )
        repository << "/" << repositoryVariantPath;
    repository << "/";
    Environment::changeRepository( _directory=boost::format( "%1%" ) % repository.str() );

    std::ostringstream meshName;
    meshName << meshCfg.name;
    auto const meshVariantTag = normalizePathFragment( options.meshVariantTag );
    if ( !meshVariantTag.empty() )
        meshName << "-" << meshVariantTag;
    meshName << "-nx" << meshCfg.nx
             << "-ny" << meshCfg.ny
             << "-nz" << meshCfg.nz
             << "-nr" << meshCfg.nr
             << "-nt" << meshCfg.nt;

    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _desc=geo( _filename=meshName.str() + ".geo",
                                           _desc=makeGeoDescription( meshCfg ),
                                           _dim=3,
                                           _order=1,
                                           _h=requestedH ? *requestedH : 1.0 ),
                                _force_rebuild=boption( "gmsh.rebuild" ) );

    return MeshBuildResult{ mesh, meshCfg.probe };
}

bool
reportFieldReferenceComparison( BenchmarkConfig::FieldReference const& reference,
                                std::string const& target,
                                shell_voigt const& epsilon,
                                shell_voigt const& sigma,
                                bool checkReference,
                                std::ostream& out )
{
    auto const it = reference.targets.find( target );
    if ( it == reference.targets.end() )
        return true;

    auto const& expected = it->second;
    double const epsilonAbsoluteError = ( epsilon - expected.epsilon ).cwiseAbs().maxCoeff();
    double const sigmaAbsoluteError = ( sigma - expected.sigma ).cwiseAbs().maxCoeff();
    double const absoluteError = std::max( epsilonAbsoluteError, sigmaAbsoluteError );
    double const referenceScale = std::max( expected.epsilon.cwiseAbs().maxCoeff(),
                                           expected.sigma.cwiseAbs().maxCoeff() );
    double const relativeError = absoluteError / std::max( referenceScale, std::numeric_limits<double>::epsilon() );
    bool const hasAbsoluteCriteria = reference.hasAbsoluteTolerance ||
                                     reference.hasEpsilonAbsoluteTolerances ||
                                     reference.hasSigmaAbsoluteTolerances;
    auto const epsilonAbsoluteTolerances =
        reference.hasEpsilonAbsoluteTolerances
            ? reference.epsilonAbsoluteTolerances
            : shell_voigt::Constant( reference.hasAbsoluteTolerance
                                         ? reference.absoluteTolerance
                                         : std::numeric_limits<double>::infinity() );
    auto const sigmaAbsoluteTolerances =
        reference.hasSigmaAbsoluteTolerances
            ? reference.sigmaAbsoluteTolerances
            : shell_voigt::Constant( reference.hasAbsoluteTolerance
                                         ? reference.absoluteTolerance
                                         : std::numeric_limits<double>::infinity() );
    bool const absoluteOk =
        ( ( epsilon - expected.epsilon ).cwiseAbs().array() <= epsilonAbsoluteTolerances.array() ).all() &&
        ( ( sigma - expected.sigma ).cwiseAbs().array() <= sigmaAbsoluteTolerances.array() ).all();

    bool ok = true;
    if ( checkReference && hasAbsoluteCriteria )
        ok = ok && absoluteOk;
    if ( checkReference && reference.hasRelativeTolerance )
        ok = ok && relativeError <= reference.relativeTolerance;

    out << "field reference '" << reference.name << "' target '" << target << "'"
        << ": abs-error=" << absoluteError
        << ", rel-error=" << relativeError;
    if ( reference.hasAbsoluteTolerance )
        out << ", abs-tol=" << reference.absoluteTolerance;
    if ( reference.hasEpsilonAbsoluteTolerances || reference.hasSigmaAbsoluteTolerances )
        out << ", component-abs-tol=yes";
    if ( reference.hasRelativeTolerance )
        out << ", rel-tol=" << reference.relativeTolerance;
    out << ", check=" << ( ok ? "PASS" : "FAIL" ) << "\n";

    bool const exceedsReportedTolerance =
        ( hasAbsoluteCriteria && !absoluteOk ) ||
        ( reference.hasRelativeTolerance && relativeError > reference.relativeTolerance );
    if ( exceedsReportedTolerance )
    {
        auto printVoigt = [&out]( char const* label, shell_voigt const& value )
        {
            out << "  " << label << " = [";
            for ( int i = 0; i < value.size(); ++i )
            {
                if ( i )
                    out << ", ";
                out << value( i );
            }
            out << "]\n";
        };
        printVoigt( "epsilon actual", epsilon );
        printVoigt( "epsilon expected", expected.epsilon );
        printVoigt( "sigma actual", sigma );
        printVoigt( "sigma expected", expected.sigma );
    }

    return ok;
}
} // namespace Feel::Quickstart::ShellBenchmark
