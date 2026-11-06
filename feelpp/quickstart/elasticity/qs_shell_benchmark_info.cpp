/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#include <feel/feelcore/environment.hpp>

#include "qs_shell_benchmark_framework.hpp"

#include <algorithm>
#include <array>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Feel;

namespace
{
namespace qsb = Feel::Quickstart::ShellBenchmark;
using BenchmarkConfig = qsb::BenchmarkConfig;
using json = qsb::json;

std::string
geometryKindName( qsb::GeometryKind kind )
{
    switch ( kind )
    {
    case qsb::GeometryKind::Rectangular:
        return "rectangular";
    case qsb::GeometryKind::Circular:
        return "circular";
    }
    return "unknown";
}

std::string
joinStrings( std::vector<std::string> const& values )
{
    if ( values.empty() )
        return "none";

    std::ostringstream os;
    for ( std::size_t i = 0; i < values.size(); ++i )
    {
        if ( i > 0 )
            os << ", ";
        os << values[i];
    }
    return os.str();
}

std::string
componentMaskString( std::array<bool, 3> const& components )
{
    static constexpr std::array<char, 3> names = { 'x', 'y', 'z' };
    std::string result;
    for ( int i = 0; i < 3; ++i )
    {
        if ( components[i] )
            result.push_back( names[i] );
    }
    return result.empty() ? "none" : result;
}

void
printIndentedLine( std::string const& label, std::string const& value, int indent = 2 )
{
    std::cout << std::string( indent, ' ' ) << label << ": " << value << "\n";
}

void
printLoads( BenchmarkConfig const& cfg )
{
    std::cout << "Loads\n";
    printIndentedLine( "body force", cfg.bodyForceExpression );

    if ( cfg.pressureLoads.empty() )
        printIndentedLine( "pressures", "none" );
    else
    {
        std::cout << "  pressures:\n";
        for ( auto const& load : cfg.pressureLoads )
            printIndentedLine( load.marker, std::to_string( load.value ), 4 );
    }

    if ( cfg.tractionLoads.empty() )
        printIndentedLine( "tractions", "none" );
    else
    {
        std::cout << "  tractions:\n";
        for ( auto const& load : cfg.tractionLoads )
            printIndentedLine( load.marker, load.expression, 4 );
    }

    if ( cfg.totalForceLoads.empty() )
        printIndentedLine( "total forces", "none" );
    else
    {
        std::cout << "  total forces:\n";
        for ( auto const& load : cfg.totalForceLoads )
            printIndentedLine( load.marker, qsb::shellVectorExpression( load.value ), 4 );
    }

    if ( cfg.pointLoads.empty() )
        printIndentedLine( "point loads", "none" );
    else
    {
        std::cout << "  point loads:\n";
        for ( auto const& load : cfg.pointLoads )
            printIndentedLine( load.marker, qsb::shellVectorExpression( load.value ), 4 );
    }
}

void
printBoundaryConditions( BenchmarkConfig const& cfg )
{
    std::cout << "Boundary Conditions\n";
    printIndentedLine( "clamp markers", joinStrings( cfg.clampMarkers ) );

    if ( cfg.pointConstraints.empty() )
    {
        printIndentedLine( "point constraints", "none" );
        return;
    }

    std::cout << "  point constraints:\n";
    for ( auto const& constraint : cfg.pointConstraints )
    {
        std::ostringstream os;
        os << "components=" << componentMaskString( constraint.components )
           << ", value=" << qsb::shellVectorExpression( constraint.value );
        printIndentedLine( constraint.marker, os.str(), 4 );
    }
}

void
printReference( BenchmarkConfig const& cfg )
{
    std::cout << "Reference\n";
    if ( !cfg.reference.hasProbeValue )
    {
        printIndentedLine( "status", "none" );
        return;
    }

    printIndentedLine( "probe value", std::to_string( cfg.reference.probeValue ) );
    if ( cfg.reference.hasAbsoluteTolerance )
        printIndentedLine( "absolute tolerance", std::to_string( cfg.reference.absoluteTolerance ) );
    if ( cfg.reference.hasRelativeTolerance )
        printIndentedLine( "relative tolerance", std::to_string( cfg.reference.relativeTolerance ) );
    if ( !cfg.reference.description.empty() )
        printIndentedLine( "description", cfg.reference.description );
}

void
printBenchmarkConfig( BenchmarkConfig const& cfg )
{
    std::cout << "Benchmark: " << cfg.name << "\n";
    if ( !cfg.description.empty() )
        std::cout << "Description: " << cfg.description << "\n";

    std::cout << "Geometry\n";
    printIndentedLine( "kind", geometryKindName( cfg.geometryKind ) );
    printIndentedLine( "template", cfg.meshTemplatePath );
    if ( cfg.geometryKind == qsb::GeometryKind::Rectangular )
    {
        printIndentedLine( "length", std::to_string( cfg.length ) );
        printIndentedLine( "width", std::to_string( cfg.width ) );
        printIndentedLine( "cells", "nx=" + std::to_string( cfg.nx ) +
                                    ", ny=" + std::to_string( cfg.ny ) +
                                    ", nz=" + std::to_string( cfg.nz ) );
    }
    else
    {
        printIndentedLine( "radius", std::to_string( cfg.radius ) );
        printIndentedLine( "cells", "nr=" + std::to_string( cfg.nr ) +
                                    ", nt=" + std::to_string( cfg.nt ) +
                                    ", nz=" + std::to_string( cfg.nz ) );
    }
    printIndentedLine( "thickness", std::to_string( cfg.thickness ) );

    std::cout << "Material\n";
    printIndentedLine( "E", std::to_string( cfg.E ) );
    printIndentedLine( "nu", std::to_string( cfg.nu ) );

    std::cout << "Probe\n";
    printIndentedLine( "label", cfg.probeLabel );
    printIndentedLine( "point", qsb::shellVectorExpression( cfg.probe ) );
    printIndentedLine( "direction", qsb::shellVectorExpression( cfg.probeDirection ) );

    printBoundaryConditions( cfg );
    printLoads( cfg );
    printReference( cfg );
}

struct BenchmarkEntry
{
    std::string name;
    std::string category;
};

void
appendBenchmarkEntries( json const& specs, char const* key, char const* category,
                        std::vector<BenchmarkEntry>& entries )
{
    if ( !specs.contains( key ) )
        return;

    for ( auto const& [name, _] : specs.at( key ).items() )
        entries.push_back( { name, category } );
}

std::vector<BenchmarkEntry>
benchmarkEntries( json const& specs )
{
    std::vector<BenchmarkEntry> entries;
    appendBenchmarkEntries( specs, "benchmarks", "benchmark", entries );
    appendBenchmarkEntries( specs, "validations", "validation", entries );

    std::sort( entries.begin(), entries.end(),
               []( BenchmarkEntry const& lhs, BenchmarkEntry const& rhs ) {
                   return std::tie( lhs.category, lhs.name ) < std::tie( rhs.category, rhs.name );
               } );
    return entries;
}

void
printBenchmarkList( json const& specs )
{
    auto const entries = benchmarkEntries( specs );
    std::cout << "Available shell cases (" << entries.size() << ")\n";
    for ( auto const& entry : entries )
    {
        auto const& spec = qsb::benchmarkSpec( specs, entry.name );
        std::cout << "  - [" << entry.category << "] " << entry.name;
        if ( spec.contains( "description" ) )
            std::cout << ": " << spec.at( "description" ).get<std::string>();
        std::cout << "\n";
    }
}

po::options_description
makeOptions()
{
    po::options_description options( "qs_shell_benchmark_info options" );
    options.add_options()
        ( "specs", po::value<std::string>(),
          "json benchmark specification file (typically provided by qs_shell_benchmark_info.cfg)" )
        ( "benchmark", po::value<std::string>()->default_value( "square-plate" ),
          "benchmark name in the JSON specification file" )
        ( "list", po::bool_switch()->default_value( false ),
          "list benchmark/validation cases available in the specification file" )
        ( "all", po::bool_switch()->default_value( false ),
          "print the full resolved configuration for every benchmark/validation case" );
    return options.add( feel_options() );
}

AboutData
makeAbout()
{
    AboutData about( "qs_shell_benchmark_info",
                     "qs_shell_benchmark_info",
                     "0.1",
                     "Display resolved shell benchmark configurations",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) Feel++ Consortium" );
    about.addAuthor( "Feel++ Consortium", "developer", "feelpp-devel@feelpp.org", "" );
    return about;
}
} // namespace

int
main( int argc, char** argv )
{
    try
    {
        Environment env( _argc=argc, _argv=argv, _desc=makeOptions(), _about=makeAbout() );

        if ( !Environment::vm().count( "specs" ) || soption( "specs" ).empty() )
            throw std::invalid_argument( "missing required option 'specs'; set it in qs_shell_benchmark_info.cfg or pass --specs" );

        auto const specsPath = Environment::expand( soption( "specs" ) );
        auto const specs = qsb::loadBenchmarkSpecs( specsPath );

        std::cout << "Specs: " << specsPath << "\n";

        if ( boption( "list" ) )
        {
            printBenchmarkList( specs );
            return 0;
        }

        if ( boption( "all" ) )
        {
            auto const entries = benchmarkEntries( specs );
            for ( std::size_t i = 0; i < entries.size(); ++i )
            {
                if ( i > 0 )
                    std::cout << "\n" << std::string( 80, '-' ) << "\n\n";
                std::cout << "Category: " << entries[i].category << "\n";
                printBenchmarkConfig( qsb::benchmarkPreset( specs, entries[i].name ) );
            }
            return 0;
        }

        printBenchmarkConfig( qsb::benchmarkPreset( specs, soption( "benchmark" ) ) );
    }
    catch ( ... )
    {
        handleExceptions();
        return 1;
    }
    return 0;
}
