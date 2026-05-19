/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_QUICKSTART_QS_ELASTICITY_CASE_HPP
#define FEELPP_QUICKSTART_QS_ELASTICITY_CASE_HPP 1

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/json.hpp>

#include "qs_elasticity_checks.hpp"

#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Feel::Quickstart::ElasticityCase
{
namespace qsec = Feel::Quickstart::ElasticityChecks;

struct MarkerExpression
{
    std::set<std::string> markers;
    std::string expression;
};

struct MarkerScalarExpression
{
    std::string marker;
    std::string expression;
};

struct PointLoadConfig
{
    std::set<std::string> markers;
    std::string expression;
    std::string quantity = "per-point";
};

struct PointConstraintConfig
{
    std::string marker;
    std::array<bool, FEELPP_DIM> components;
    std::array<std::string, FEELPP_DIM> values;
};

struct CheckerConfig
{
    std::optional<bool> check;
    std::optional<bool> exact;
    std::optional<std::string> solution;
    std::optional<std::string> script;
    std::optional<double> exactTolerance;
    std::optional<double> orderTolerance;
};

struct Config
{
    double E = 1.0e6;
    double nu = 0.3;
    std::string meshFilename;
    std::optional<size_type> expectedElements;
    std::string bodyForceExpression;
    std::string dirichletExpression;
    std::set<std::string> dirichletMarkers{ "Dirichlet" };
    std::vector<MarkerExpression> faceTractions;
    std::vector<MarkerExpression> faceTotalForces;
    std::vector<MarkerScalarExpression> facePressures;
    std::vector<PointLoadConfig> pointForces;
    std::vector<PointLoadConfig> pointMoments;
    std::vector<PointConstraintConfig> pointConstraints;
    CheckerConfig checker;
    qsec::Config<FEELPP_DIM> referenceChecks;
};

inline std::string
zeroVectorExpression()
{
    if constexpr ( FEELPP_DIM == 2 )
        return "{0,0}";
    else
        return "{0,0,0}";
}

inline std::string
zeroMomentExpression()
{
    if constexpr ( FEELPP_DIM == 2 )
        return "0";
    else
        return "{0,0,0}";
}

inline Config
defaultConfig( std::string const& zeroVectorExpr = zeroVectorExpression() )
{
    Config cfg;
    cfg.bodyForceExpression = zeroVectorExpr;
    cfg.dirichletExpression = zeroVectorExpr;
    return cfg;
}

inline void
addOptions( po::options_description& options,
            std::string const& zeroVectorExpr = zeroVectorExpression(),
            std::string const& zeroMomentExpr = zeroMomentExpression() )
{
    options.add_options()
        ( "E", po::value<double>()->default_value( 1.0e6 ), "Young modulus" )
        ( "nu", po::value<double>()->default_value( 0.3 ), "Poisson ratio" )
        ( "elasticity.case", po::value<std::string>()->default_value( "" ), "single JSON elasticity case file" )
        ( "checks.target", po::value<std::string>()->default_value( "" ), "reference target used by JSON checks" )
        ( "mesh.expected-elements", po::value<int>()->default_value( -1 ), "fail if the loaded mesh does not have this global number of elements" )
        ( "no-solve", po::value<bool>()->default_value( false ), "No solve" )
        ( "weakdir", po::value<bool>()->default_value( false ), "use weak dirichlet" )
        ( "gamma", po::value<double>()->default_value( 100 ), "penalisation term" )
        ( "nullspace", po::value<bool>()->default_value( false ), "add null space" )
        ( "dirichlet.markers", po::value<std::vector<std::string>>()->multitoken(), "face markers receiving the displacement Dirichlet condition" )
        ( "point-force.markers", po::value<std::vector<std::string>>()->multitoken(), "point markers receiving a point force" )
        ( "point-force.expr", po::value<std::string>()->default_value( zeroVectorExpr ), "point force vector expression" )
        ( "point-force.quantity", po::value<std::string>()->default_value( "per-point" ), "point force quantity: per-point or total" )
        ( "point-moment.markers", po::value<std::vector<std::string>>()->multitoken(), "point markers receiving a point moment" )
        ( "point-moment.expr", po::value<std::string>()->default_value( zeroMomentExpr ), "point moment expression: scalar in 2D, vector in 3D" )
        ( "point-moment.quantity", po::value<std::string>()->default_value( "per-point" ), "point moment quantity: per-point or total" )
        ( "cantilever.check", po::value<bool>()->default_value( false ), "check tip displacement against an Euler-Bernoulli cantilever reference" )
        ( "cantilever.tip-marker", po::value<std::string>()->default_value( "tip" ), "point marker used for the cantilever tip displacement check" )
        ( "cantilever.component", po::value<int>()->default_value( FEELPP_DIM-1 ), "displacement component used by the cantilever check" )
        ( "cantilever.length", po::value<double>()->default_value( 1. ), "cantilever beam length" )
        ( "cantilever.height", po::value<double>()->default_value( 1. ), "cantilever bending height" )
        ( "cantilever.thickness", po::value<double>()->default_value( 1. ), "cantilever out-of-plane thickness" )
        ( "cantilever.tip-force", po::value<double>()->default_value( 0. ), "signed transverse tip force used in the cantilever reference" )
        ( "cantilever.tip-moment", po::value<double>()->default_value( 0. ), "signed end moment used in the cantilever reference" )
        ( "cantilever.tolerance.relative", po::value<double>()->default_value( 0.5 ), "relative tolerance for the cantilever check" )
        ( "cantilever.tolerance.absolute", po::value<double>()->default_value( 1e-12 ), "absolute tolerance for the cantilever check" )
        ;
}

inline bool
optionExplicitlySet( std::string const& name )
{
    auto const& vm = Feel::Environment::vm();
    auto const it = vm.find( name );
    return it != vm.end() && !it->second.defaulted();
}

inline std::string
formatDouble( double value )
{
    std::ostringstream os;
    os << std::setprecision( 16 ) << value;
    return os.str();
}

inline std::string
coordinateSuffix()
{
    if constexpr ( FEELPP_DIM == 2 )
        return ":x:y";
    else
        return ":x:y:z";
}

inline std::string
jsonExpression( Feel::nl::json const& value, std::string const& description )
{
    if ( value.is_string() )
        return value.get<std::string>();
    if ( value.is_number() )
        return formatDouble( value.get<double>() );
    if ( value.is_array() )
    {
        std::ostringstream os;
        os << "{";
        for ( std::size_t i = 0; i < value.size(); ++i )
        {
            if ( i > 0 )
                os << ",";
            auto const& component = value.at( i );
            if ( component.is_string() )
                os << component.get<std::string>();
            else if ( component.is_number() )
                os << formatDouble( component.get<double>() );
            else
                throw std::invalid_argument( description + " contains a non-scalar component" );
        }
        os << "}" << coordinateSuffix();
        return os.str();
    }

    throw std::invalid_argument( description + " must be a scalar, string expression, or array expression" );
}

inline std::string
jsonExpression( Feel::nl::json const& object,
                std::string const& key,
                std::string const& defaultValue )
{
    if ( !object.contains( key ) )
        return defaultValue;
    return jsonExpression( object.at( key ), key );
}

inline Feel::nl::json const&
jsonObject( Feel::nl::json const& object, std::string const& key )
{
    if ( !object.contains( key ) || !object.at( key ).is_object() )
        throw std::invalid_argument( "missing json object '" + key + "'" );
    return object.at( key );
}

inline std::filesystem::path
caseRelativePath( std::string const& value, std::filesystem::path const& caseDirectory )
{
    auto expanded = Feel::Environment::expand( value );
    std::filesystem::path path( expanded );
    if ( path.is_absolute() || expanded.find( '$' ) != std::string::npos )
        return path.lexically_normal();
    return ( caseDirectory / path ).lexically_normal();
}

inline std::set<std::string>
jsonMarkers( Feel::nl::json const& object, std::string const& context )
{
    std::set<std::string> markers;
    if ( object.contains( "marker" ) )
        markers.insert( object.at( "marker" ).get<std::string>() );
    if ( object.contains( "markers" ) )
    {
        for ( auto const& marker : object.at( "markers" ) )
            markers.insert( marker.get<std::string>() );
    }
    if ( markers.empty() )
        throw std::invalid_argument( context + " must define 'marker' or 'markers'" );
    return markers;
}

inline std::array<bool, FEELPP_DIM>
jsonComponentMask( Feel::nl::json const& object, std::string const& key, bool defaultValue )
{
    std::array<bool, FEELPP_DIM> components;
    components.fill( defaultValue );
    if ( !object.contains( key ) )
        return components;

    auto const& values = object.at( key );
    if ( !values.is_array() || values.size() < FEELPP_DIM )
        throw std::invalid_argument( key + " must be a boolean array with at least " + std::to_string( FEELPP_DIM ) + " entries" );
    for ( int c = 0; c < FEELPP_DIM; ++c )
        components[c] = values.at( c ).get<bool>();
    return components;
}

inline std::array<std::string, FEELPP_DIM>
jsonComponentExpressions( Feel::nl::json const& object, std::string const& key )
{
    std::array<std::string, FEELPP_DIM> values;
    values.fill( "0" );
    if ( !object.contains( key ) )
        return values;

    auto const& value = object.at( key );
    if ( value.is_array() )
    {
        if ( value.size() < FEELPP_DIM )
            throw std::invalid_argument( key + " must contain at least " + std::to_string( FEELPP_DIM ) + " components" );
        for ( int c = 0; c < FEELPP_DIM; ++c )
            values[c] = jsonExpression( value.at( c ), key );
        return values;
    }

    auto const scalarValue = jsonExpression( value, key );
    values.fill( scalarValue );
    return values;
}

inline Feel::nl::json const&
singleCaseObject( Feel::nl::json const& root )
{
    if ( root.contains( "case" ) )
        return jsonObject( root, "case" );

    for ( auto const* key : { "benchmarks", "validations" } )
    {
        if ( root.contains( key ) )
        {
            auto const& cases = jsonObject( root, key );
            if ( cases.size() != 1 )
                throw std::invalid_argument( std::string( "json '" ) + key + "' catalog must be split into one file per qs_elasticity test" );
            return cases.begin().value();
        }
    }

    return root;
}

inline void
applyJson( Config& cfg,
           std::string const& casePath,
           std::string const& zeroVectorExpression,
           std::string const& zeroMomentExpression )
{
    auto const resolvedPath = std::filesystem::path( Feel::Environment::expand( casePath ) ).lexically_normal();
    std::ifstream input( resolvedPath );
    if ( !input )
        throw std::invalid_argument( "unable to open elasticity case json '" + resolvedPath.string() + "'" );

    auto const root = Feel::nl::json::parse( input );
    auto const& specs = singleCaseObject( root );
    auto const caseDirectory = resolvedPath.parent_path();
    cfg.referenceChecks = qsec::jsonConfig<FEELPP_DIM>( specs );

    if ( specs.contains( "material" ) )
    {
        auto const& material = jsonObject( specs, "material" );
        cfg.E = material.value( "E", cfg.E );
        cfg.nu = material.value( "nu", cfg.nu );
    }

    if ( specs.contains( "gmsh" ) && specs.at( "gmsh" ).contains( "filename" ) )
        cfg.meshFilename = caseRelativePath( specs.at( "gmsh" ).at( "filename" ).get<std::string>(), caseDirectory ).string();
    if ( specs.contains( "mesh" ) && specs.at( "mesh" ).contains( "filename" ) )
        cfg.meshFilename = caseRelativePath( specs.at( "mesh" ).at( "filename" ).get<std::string>(), caseDirectory ).string();
    if ( specs.contains( "mesh" ) && specs.at( "mesh" ).contains( "expectedElements" ) )
        cfg.expectedElements = specs.at( "mesh" ).at( "expectedElements" ).get<size_type>();

    if ( specs.contains( "functions" ) )
    {
        auto const& functions = jsonObject( specs, "functions" );
        cfg.bodyForceExpression = jsonExpression( functions, "f", cfg.bodyForceExpression );
        cfg.dirichletExpression = jsonExpression( functions, "g", cfg.dirichletExpression );
    }

    if ( specs.contains( "boundaryConditions" ) )
    {
        auto const& boundary = jsonObject( specs, "boundaryConditions" );
        if ( boundary.contains( "clampMarkers" ) )
            cfg.dirichletMarkers = boundary.at( "clampMarkers" ).get<std::set<std::string>>();
        if ( boundary.contains( "dirichletMarkers" ) )
            cfg.dirichletMarkers = boundary.at( "dirichletMarkers" ).get<std::set<std::string>>();
        if ( boundary.contains( "dirichlet" ) )
        {
            auto const& dirichlet = jsonObject( boundary, "dirichlet" );
            cfg.dirichletMarkers = jsonMarkers( dirichlet, "boundaryConditions.dirichlet" );
            cfg.dirichletExpression = jsonExpression( dirichlet, "value", cfg.dirichletExpression );
        }
        if ( boundary.contains( "pointConstraints" ) )
        {
            for ( auto const& constraint : boundary.at( "pointConstraints" ) )
            {
                cfg.pointConstraints.push_back(
                    { constraint.at( "marker" ).get<std::string>(),
                      jsonComponentMask( constraint, "components", true ),
                      jsonComponentExpressions( constraint, "value" ) } );
            }
        }
    }

    if ( specs.contains( "loads" ) )
    {
        auto const& loads = jsonObject( specs, "loads" );
        cfg.bodyForceExpression = jsonExpression( loads, "bodyForce", cfg.bodyForceExpression );

        if ( loads.contains( "tractions" ) )
        {
            for ( auto const& load : loads.at( "tractions" ) )
                cfg.faceTractions.push_back( { jsonMarkers( load, "loads.tractions[]" ),
                                               jsonExpression( load, "expr", jsonExpression( load, "value", zeroVectorExpression ) ) } );
        }

        if ( loads.contains( "pressures" ) )
        {
            for ( auto const& load : loads.at( "pressures" ) )
            {
                auto const markers = jsonMarkers( load, "loads.pressures[]" );
                auto const pressure = jsonExpression( load, "value", "0" );
                for ( auto const& marker : markers )
                    cfg.facePressures.push_back( { marker, pressure } );
            }
        }

        if ( loads.contains( "totalForces" ) )
        {
            for ( auto const& load : loads.at( "totalForces" ) )
                cfg.faceTotalForces.push_back( { jsonMarkers( load, "loads.totalForces[]" ),
                                                 jsonExpression( load, "value", zeroVectorExpression ) } );
        }

        if ( loads.contains( "pointLoads" ) )
        {
            for ( auto const& load : loads.at( "pointLoads" ) )
                cfg.pointForces.push_back( { jsonMarkers( load, "loads.pointLoads[]" ),
                                             jsonExpression( load, "value", zeroVectorExpression ),
                                             load.value( "quantity", "per-point" ) } );
        }
        if ( loads.contains( "pointForces" ) )
        {
            for ( auto const& load : loads.at( "pointForces" ) )
                cfg.pointForces.push_back( { jsonMarkers( load, "loads.pointForces[]" ),
                                             jsonExpression( load, "value", zeroVectorExpression ),
                                             load.value( "quantity", "per-point" ) } );
        }
        if ( loads.contains( "pointMoments" ) )
        {
            for ( auto const& load : loads.at( "pointMoments" ) )
                cfg.pointMoments.push_back( { jsonMarkers( load, "loads.pointMoments[]" ),
                                              jsonExpression( load, "value", zeroMomentExpression ),
                                              load.value( "quantity", "per-point" ) } );
        }
    }

    if ( specs.contains( "checker" ) )
    {
        auto const& checker = jsonObject( specs, "checker" );
        if ( checker.contains( "check" ) )
            cfg.checker.check = checker.at( "check" ).get<bool>();
        if ( checker.contains( "exact" ) )
            cfg.checker.exact = checker.at( "exact" ).get<bool>();
        if ( checker.contains( "solution" ) )
            cfg.checker.solution = jsonExpression( checker.at( "solution" ), "checker.solution" );
        if ( checker.contains( "script" ) )
            cfg.checker.script = caseRelativePath( checker.at( "script" ).get<std::string>(), caseDirectory ).string();
        if ( checker.contains( "tolerance" ) )
        {
            auto const& tolerance = jsonObject( checker, "tolerance" );
            if ( tolerance.contains( "exact" ) )
                cfg.checker.exactTolerance = tolerance.at( "exact" ).get<double>();
            if ( tolerance.contains( "order" ) )
                cfg.checker.orderTolerance = tolerance.at( "order" ).get<double>();
        }
        if ( checker.contains( "tolerance.exact" ) )
            cfg.checker.exactTolerance = checker.at( "tolerance.exact" ).get<double>();
        if ( checker.contains( "tolerance.order" ) )
            cfg.checker.orderTolerance = checker.at( "tolerance.order" ).get<double>();
    }
}

inline std::set<std::string>
markersFromOption( std::string const& opt )
{
    std::set<std::string> markers;
    if ( Feel::Environment::vm().count( opt ) )
    {
        auto markerList = Feel::Environment::vm()[opt].as<std::vector<std::string>>();
        markers.insert( markerList.begin(), markerList.end() );
    }
    return markers;
}

inline void
applyEnvironmentOverrides( Config& cfg )
{
    if ( optionExplicitlySet( "E" ) )
        cfg.E = doption( "E" );
    if ( optionExplicitlySet( "nu" ) )
        cfg.nu = doption( "nu" );
    if ( optionExplicitlySet( "gmsh.filename" ) )
        cfg.meshFilename = soption( "gmsh.filename" );
    if ( optionExplicitlySet( "mesh.expected-elements" ) && ioption( "mesh.expected-elements" ) >= 0 )
        cfg.expectedElements = static_cast<size_type>( ioption( "mesh.expected-elements" ) );
    if ( optionExplicitlySet( "checks.target" ) && !soption( "checks.target" ).empty() )
        cfg.referenceChecks.target = soption( "checks.target" );
    if ( optionExplicitlySet( "functions.f" ) )
        cfg.bodyForceExpression = soption( "functions.f" );
    if ( optionExplicitlySet( "functions.g" ) )
        cfg.dirichletExpression = soption( "functions.g" );
    if ( auto markers = markersFromOption( "dirichlet.markers" ); !markers.empty() )
        cfg.dirichletMarkers = std::move( markers );
    if ( auto markers = markersFromOption( "point-force.markers" ); !markers.empty() )
        cfg.pointForces.push_back( { std::move( markers ), soption( "point-force.expr" ), soption( "point-force.quantity" ) } );
    if ( auto markers = markersFromOption( "point-moment.markers" ); !markers.empty() )
        cfg.pointMoments.push_back( { std::move( markers ), soption( "point-moment.expr" ), soption( "point-moment.quantity" ) } );
    if ( optionExplicitlySet( "checker.check" ) )
        cfg.checker.check = boption( "checker.check" );
    if ( optionExplicitlySet( "checker.exact" ) )
        cfg.checker.exact = boption( "checker.exact" );
    if ( optionExplicitlySet( "checker.solution" ) )
        cfg.checker.solution = soption( "checker.solution" );
    if ( optionExplicitlySet( "checker.script" ) )
        cfg.checker.script = soption( "checker.script" );
    if ( optionExplicitlySet( "checker.tolerance.exact" ) )
        cfg.checker.exactTolerance = doption( "checker.tolerance.exact" );
    if ( optionExplicitlySet( "checker.tolerance.order" ) )
        cfg.checker.orderTolerance = doption( "checker.tolerance.order" );

    if ( boption( "cantilever.check" ) )
    {
        qsec::CantileverReference<FEELPP_DIM> cantilever;
        cantilever.name = "command-line cantilever";
        cantilever.hasMarker = true;
        cantilever.tipMarker = soption( "cantilever.tip-marker" );
        cantilever.component = ioption( "cantilever.component" );
        cantilever.length = doption( "cantilever.length" );
        cantilever.height = doption( "cantilever.height" );
        cantilever.thickness = doption( "cantilever.thickness" );
        cantilever.tipForce = doption( "cantilever.tip-force" );
        cantilever.tipMoment = doption( "cantilever.tip-moment" );
        cantilever.tolerance.hasRelative = true;
        cantilever.tolerance.relative = doption( "cantilever.tolerance.relative" );
        cantilever.tolerance.hasAbsolute = true;
        cantilever.tolerance.absolute = doption( "cantilever.tolerance.absolute" );
        cfg.referenceChecks.cantilevers.push_back( std::move( cantilever ) );
    }
}

inline Config
fromEnvironment( std::string const& zeroVectorExpr = zeroVectorExpression(),
                 std::string const& zeroMomentExpr = zeroMomentExpression() )
{
    auto cfg = defaultConfig( zeroVectorExpr );
    if ( optionExplicitlySet( "elasticity.case" ) && !soption( "elasticity.case" ).empty() )
        applyJson( cfg, soption( "elasticity.case" ), zeroVectorExpr, zeroMomentExpr );
    applyEnvironmentOverrides( cfg );
    return cfg;
}

} // namespace Feel::Quickstart::ElasticityCase

#endif
