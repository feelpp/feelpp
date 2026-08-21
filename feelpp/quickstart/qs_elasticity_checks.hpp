/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_QUICKSTART_QS_ELASTICITY_CHECKS_HPP
#define FEELPP_QUICKSTART_QS_ELASTICITY_CHECKS_HPP 1

#include <feel/feelcore/checker.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feeldiscr/pdhm.hpp>
#include <feel/feelvf/vf.hpp>

#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Feel::Quickstart::ElasticityChecks
{
template<int Dim>
using point_type = Eigen::Matrix<double, Dim, 1>;

using voigt_type = Eigen::Matrix<double, 6, 1>;

struct Tolerance
{
    bool hasAbsolute = false;
    double absolute = 0.0;
    bool hasRelative = false;
    double relative = 0.0;
};

inline bool
hasTolerance( Tolerance const& tolerance )
{
    return tolerance.hasAbsolute || tolerance.hasRelative;
}

struct TargetCheckOptions
{
    bool hasCheck = false;
    bool check = true;
    Tolerance tolerance;
};

template<int Dim>
struct ProbeReference
{
    bool hasProbe = false;
    point_type<Dim> point = point_type<Dim>::Zero();
    point_type<Dim> direction = point_type<Dim>::Zero();
    std::string label = "probe displacement";

    bool hasReference = false;
    bool check = true;
    double value = 0.0;
    Tolerance tolerance;
    std::string description;
};

template<int Dim>
struct CantileverReference
{
    std::string name = "Euler-Bernoulli cantilever";
    std::string description;
    bool check = true;
    Tolerance tolerance;

    bool hasPoint = false;
    point_type<Dim> point = point_type<Dim>::Zero();
    point_type<Dim> direction = point_type<Dim>::Zero();
    bool hasMarker = false;
    std::string tipMarker = "tip";
    int component = Dim - 1;

    double length = 1.0;
    double height = 1.0;
    double thickness = 1.0;
    double tipForce = 0.0;
    double tipMoment = 0.0;
    std::map<std::string, TargetCheckOptions> targets;
};

struct FieldReferenceValues
{
    bool hasEpsilon = false;
    voigt_type epsilon = voigt_type::Zero();
    bool hasSigma = false;
    voigt_type sigma = voigt_type::Zero();
};

template<int Dim>
struct FieldReference
{
    std::string name = "field-reference";
    std::string description;
    bool check = true;
    point_type<Dim> point = point_type<Dim>::Zero();
    Tolerance tolerance;
    std::map<std::string, FieldReferenceValues> targets;
};

template<int Dim>
struct Config
{
    std::string target = "hexa8";
    std::vector<ProbeReference<Dim>> probes;
    std::vector<CantileverReference<Dim>> cantilevers;
    std::vector<FieldReference<Dim>> fieldReferences;

    bool hasProbeReference() const
    {
        return std::any_of( probes.begin(), probes.end(), []( auto const& probe ) {
            return probe.hasProbe && probe.hasReference;
        } );
    }
    bool hasCantileverReferences() const { return !cantilevers.empty(); }
    bool hasFieldReferences() const { return !fieldReferences.empty(); }
    bool hasActiveChecks() const { return hasProbeReference() || hasCantileverReferences() || hasFieldReferences(); }
};

inline std::string
formatDouble( double value )
{
    std::ostringstream os;
    os << std::setprecision( 16 ) << value;
    return os.str();
}

inline double
relativeError( double computed, double reference )
{
    double const absoluteError = computed - reference;
    return std::abs( reference ) > 0.0 ? std::abs( absoluteError/reference ) : std::abs( absoluteError );
}

inline double
jsonDouble( nl::json const& value, std::string const& context )
{
    if ( !value.is_number() )
        throw std::invalid_argument( context + " must be numeric" );
    return value.get<double>();
}

template<int Dim>
point_type<Dim>
jsonPoint( nl::json const& value, std::string const& context )
{
    if ( !value.is_array() || value.size() < Dim )
        throw std::invalid_argument( context + " must be a numeric array with at least " + std::to_string( Dim ) + " entries" );

    point_type<Dim> point = point_type<Dim>::Zero();
    for ( int i = 0; i < Dim; ++i )
        point( i ) = jsonDouble( value.at( i ), context );
    return point;
}

inline voigt_type
jsonVoigt( nl::json const& value, std::string const& context )
{
    if ( !value.is_array() || value.size() < 6 )
        throw std::invalid_argument( context + " must be a numeric array with at least six entries" );

    voigt_type result = voigt_type::Zero();
    for ( int i = 0; i < 6; ++i )
        result( i ) = jsonDouble( value.at( i ), context );
    return result;
}

inline Tolerance
jsonTolerance( nl::json const& object )
{
    Tolerance tolerance;
    if ( object.contains( "tolerance" ) )
    {
        auto const& nested = object.at( "tolerance" );
        if ( !nested.is_object() )
            throw std::invalid_argument( "tolerance must be an object" );
        if ( nested.contains( "absolute" ) )
        {
            tolerance.hasAbsolute = true;
            tolerance.absolute = nested.at( "absolute" ).get<double>();
        }
        if ( nested.contains( "absoluteTolerance" ) )
        {
            tolerance.hasAbsolute = true;
            tolerance.absolute = nested.at( "absoluteTolerance" ).get<double>();
        }
        if ( nested.contains( "relative" ) )
        {
            tolerance.hasRelative = true;
            tolerance.relative = nested.at( "relative" ).get<double>();
        }
        if ( nested.contains( "relativeTolerance" ) )
        {
            tolerance.hasRelative = true;
            tolerance.relative = nested.at( "relativeTolerance" ).get<double>();
        }
    }
    if ( object.contains( "absoluteTolerance" ) )
    {
        tolerance.hasAbsolute = true;
        tolerance.absolute = object.at( "absoluteTolerance" ).get<double>();
    }
    if ( object.contains( "relativeTolerance" ) )
    {
        tolerance.hasRelative = true;
        tolerance.relative = object.at( "relativeTolerance" ).get<double>();
    }
    return tolerance;
}

inline double
jsonObjectValue( nl::json const& object,
                 std::string const& key,
                 double defaultValue )
{
    return object.contains( key ) ? object.at( key ).get<double>() : defaultValue;
}

inline TargetCheckOptions
jsonTargetCheckOptions( nl::json const& object, std::string const& context )
{
    TargetCheckOptions options;
    if ( object.is_boolean() )
    {
        options.hasCheck = true;
        options.check = object.get<bool>();
        return options;
    }
    if ( !object.is_object() )
        throw std::invalid_argument( context + " must be a boolean or an object" );

    if ( object.contains( "check" ) )
    {
        options.hasCheck = true;
        options.check = object.at( "check" ).get<bool>();
    }
    options.tolerance = jsonTolerance( object );
    return options;
}

template<int Dim>
Config<Dim>
jsonConfig( nl::json const& specs )
{
    Config<Dim> config;

    auto parseCantileverReference = []( nl::json const& item, std::string const& context )
    {
        CantileverReference<Dim> cantilever;
        if ( item.contains( "name" ) )
            cantilever.name = item.at( "name" ).get<std::string>();
        if ( item.contains( "label" ) )
            cantilever.name = item.at( "label" ).get<std::string>();
        if ( item.contains( "description" ) )
            cantilever.description = item.at( "description" ).get<std::string>();
        if ( item.contains( "check" ) )
        {
            if ( item.at( "check" ).is_boolean() )
                cantilever.check = item.at( "check" ).get<bool>();
            else if ( item.at( "check" ).is_object() )
            {
                for ( auto const& [targetName, targetSpec] : item.at( "check" ).items() )
                    cantilever.targets[targetName] = jsonTargetCheckOptions( targetSpec, context + ".check." + targetName );
            }
            else
                throw std::invalid_argument( context + ".check must be a boolean or target object" );
        }
        cantilever.tolerance = jsonTolerance( item );

        if ( item.contains( "component" ) )
            cantilever.component = item.at( "component" ).get<int>();
        if ( cantilever.component < 0 || cantilever.component >= Dim )
            throw std::invalid_argument( context + ".component must be in [0," + std::to_string( Dim-1 ) + "]" );

        if ( item.contains( "point" ) )
        {
            cantilever.hasPoint = true;
            cantilever.point = jsonPoint<Dim>( item.at( "point" ), context + ".point" );
        }
        if ( item.contains( "direction" ) )
            cantilever.direction = jsonPoint<Dim>( item.at( "direction" ), context + ".direction" );
        else
            cantilever.direction( cantilever.component ) = 1.0;

        if ( item.contains( "tipMarker" ) )
        {
            cantilever.hasMarker = true;
            cantilever.tipMarker = item.at( "tipMarker" ).get<std::string>();
        }
        if ( item.contains( "marker" ) )
        {
            cantilever.hasMarker = true;
            cantilever.tipMarker = item.at( "marker" ).get<std::string>();
        }

        if ( item.contains( "tip" ) )
        {
            auto const& tip = item.at( "tip" );
            if ( !tip.is_object() )
                throw std::invalid_argument( context + ".tip must be an object" );
            if ( tip.contains( "point" ) )
            {
                cantilever.hasPoint = true;
                cantilever.point = jsonPoint<Dim>( tip.at( "point" ), context + ".tip.point" );
            }
            if ( tip.contains( "direction" ) )
                cantilever.direction = jsonPoint<Dim>( tip.at( "direction" ), context + ".tip.direction" );
            if ( tip.contains( "component" ) )
                cantilever.component = tip.at( "component" ).get<int>();
            if ( cantilever.component < 0 || cantilever.component >= Dim )
                throw std::invalid_argument( context + ".tip.component must be in [0," + std::to_string( Dim-1 ) + "]" );
            if ( !tip.contains( "direction" ) )
            {
                cantilever.direction = point_type<Dim>::Zero();
                cantilever.direction( cantilever.component ) = 1.0;
            }
            if ( tip.contains( "marker" ) )
            {
                cantilever.hasMarker = true;
                cantilever.tipMarker = tip.at( "marker" ).get<std::string>();
            }
        }

        if ( item.contains( "beam" ) )
        {
            auto const& beam = item.at( "beam" );
            if ( !beam.is_object() )
                throw std::invalid_argument( context + ".beam must be an object" );
            cantilever.length = jsonObjectValue( beam, "length", cantilever.length );
            cantilever.height = jsonObjectValue( beam, "height", cantilever.height );
            cantilever.thickness = jsonObjectValue( beam, "thickness", cantilever.thickness );
        }
        cantilever.length = jsonObjectValue( item, "length", cantilever.length );
        cantilever.height = jsonObjectValue( item, "height", cantilever.height );
        cantilever.thickness = jsonObjectValue( item, "thickness", cantilever.thickness );

        if ( item.contains( "load" ) )
        {
            auto const& load = item.at( "load" );
            if ( !load.is_object() )
                throw std::invalid_argument( context + ".load must be an object" );
            cantilever.tipForce = jsonObjectValue( load, "tipForce", cantilever.tipForce );
            cantilever.tipForce = jsonObjectValue( load, "force", cantilever.tipForce );
            cantilever.tipMoment = jsonObjectValue( load, "tipMoment", cantilever.tipMoment );
            cantilever.tipMoment = jsonObjectValue( load, "moment", cantilever.tipMoment );
        }
        cantilever.tipForce = jsonObjectValue( item, "tipForce", cantilever.tipForce );
        cantilever.tipMoment = jsonObjectValue( item, "tipMoment", cantilever.tipMoment );

        if ( item.contains( "targets" ) )
        {
            if ( !item.at( "targets" ).is_object() )
                throw std::invalid_argument( context + ".targets must be an object" );
            for ( auto const& [targetName, targetSpec] : item.at( "targets" ).items() )
                cantilever.targets[targetName] = jsonTargetCheckOptions( targetSpec, context + ".targets." + targetName );
        }

        if ( !cantilever.hasPoint && !cantilever.hasMarker )
            cantilever.hasMarker = true;
        return cantilever;
    };

    auto parseFieldReference = []( nl::json const& fieldSpec, std::string const& context )
    {
        FieldReference<Dim> fieldReference;
        if ( fieldSpec.contains( "name" ) )
            fieldReference.name = fieldSpec.at( "name" ).get<std::string>();
        if ( fieldSpec.contains( "description" ) )
            fieldReference.description = fieldSpec.at( "description" ).get<std::string>();
        fieldReference.check = fieldSpec.value( "check", true );
        if ( !fieldSpec.contains( "point" ) )
            throw std::invalid_argument( context + ".point is required" );
        fieldReference.point = jsonPoint<Dim>( fieldSpec.at( "point" ), context + ".point" );
        fieldReference.tolerance = jsonTolerance( fieldSpec );

        char const* targetKey = fieldSpec.contains( "references" ) ? "references" : "targets";
        if ( !fieldSpec.contains( targetKey ) || !fieldSpec.at( targetKey ).is_object() )
            throw std::invalid_argument( context + "." + targetKey + " must be an object" );
        for ( auto const& [targetName, targetSpec] : fieldSpec.at( targetKey ).items() )
        {
            FieldReferenceValues values;
            if ( targetSpec.contains( "epsilon" ) )
            {
                values.hasEpsilon = true;
                values.epsilon = jsonVoigt( targetSpec.at( "epsilon" ), context + "." + targetKey + "." + targetName + ".epsilon" );
            }
            if ( targetSpec.contains( "sigma" ) )
            {
                values.hasSigma = true;
                values.sigma = jsonVoigt( targetSpec.at( "sigma" ), context + "." + targetKey + "." + targetName + ".sigma" );
            }
            if ( !values.hasEpsilon && !values.hasSigma )
                throw std::invalid_argument( context + "." + targetKey + "." + targetName + " must define epsilon or sigma" );
            fieldReference.targets.emplace( targetName, values );
        }
        return fieldReference;
    };

    if ( specs.contains( "checks" ) )
    {
        auto const& checks = specs.at( "checks" );
        if ( !checks.is_object() )
            throw std::invalid_argument( "checks must be an object" );
        config.target = checks.value( "target", config.target );
        if ( !checks.contains( "items" ) )
            return config;
        if ( !checks.at( "items" ).is_array() )
            throw std::invalid_argument( "checks.items must be an array" );

        for ( auto const& item : checks.at( "items" ) )
        {
            std::string const type = item.value( "type", std::string{} );
            if ( type == "probe" || type == "displacement-probe" )
            {
                ProbeReference<Dim> probe;
                probe.hasProbe = true;
                if ( item.contains( "name" ) )
                    probe.label = item.at( "name" ).get<std::string>();
                if ( item.contains( "label" ) )
                    probe.label = item.at( "label" ).get<std::string>();
                if ( item.contains( "description" ) )
                    probe.description = item.at( "description" ).get<std::string>();
                probe.check = item.value( "check", true );
                probe.point = jsonPoint<Dim>( item.at( "point" ), "checks.items[].point" );
                probe.direction = jsonPoint<Dim>( item.at( "direction" ), "checks.items[].direction" );
                probe.tolerance = jsonTolerance( item );

                if ( item.contains( "reference" ) )
                {
                    probe.hasReference = true;
                    auto const& reference = item.at( "reference" );
                    if ( reference.is_number() )
                        probe.value = reference.get<double>();
                    else if ( reference.is_object() )
                    {
                        probe.check = reference.value( "check", probe.check );
                        if ( reference.contains( "value" ) )
                            probe.value = reference.at( "value" ).get<double>();
                        else if ( reference.contains( "probeValue" ) )
                            probe.value = reference.at( "probeValue" ).get<double>();
                        else
                            throw std::invalid_argument( "checks.items[].reference must define value" );
                        probe.tolerance = jsonTolerance( reference );
                        if ( reference.contains( "description" ) )
                            probe.description = reference.at( "description" ).get<std::string>();
                    }
                    else
                        throw std::invalid_argument( "checks.items[].reference must be numeric or an object" );
                    config.probes.push_back( std::move( probe ) );
                }
            }
            else if ( type == "elasticity-point-fields" || type == "field-reference" )
            {
                config.fieldReferences.push_back( parseFieldReference( item, "checks.items[]" ) );
            }
            else if ( type == "euler-bernoulli-cantilever" || type == "cantilever-euler-bernoulli" || type == "cantilever" )
            {
                config.cantilevers.push_back( parseCantileverReference( item, "checks.items[]" ) );
            }
            else
                throw std::invalid_argument( "unsupported checks.items[] type '" + type + "'" );
        }
        return config;
    }

    ProbeReference<Dim> legacyProbe;
    if ( specs.contains( "probe" ) )
    {
        auto const& probe = specs.at( "probe" );
        if ( !probe.is_object() )
            throw std::invalid_argument( "probe must be an object" );
        if ( probe.contains( "point" ) )
        {
            legacyProbe.hasProbe = true;
            legacyProbe.point = jsonPoint<Dim>( probe.at( "point" ), "probe.point" );
        }
        if ( probe.contains( "direction" ) )
            legacyProbe.direction = jsonPoint<Dim>( probe.at( "direction" ), "probe.direction" );
        if ( probe.contains( "label" ) )
            legacyProbe.label = probe.at( "label" ).get<std::string>();
    }

    if ( specs.contains( "reference" ) )
    {
        auto const& reference = specs.at( "reference" );
        if ( !reference.is_object() )
            throw std::invalid_argument( "reference must be an object" );
        if ( reference.contains( "probeValue" ) )
        {
            if ( !legacyProbe.hasProbe )
                throw std::invalid_argument( "reference.probeValue requires probe.point" );
            legacyProbe.hasReference = true;
            legacyProbe.check = reference.value( "check", true );
            legacyProbe.value = reference.at( "probeValue" ).get<double>();
            legacyProbe.tolerance = jsonTolerance( reference );
            if ( reference.contains( "description" ) )
                legacyProbe.description = reference.at( "description" ).get<std::string>();
            config.probes.push_back( std::move( legacyProbe ) );
        }
    }

    if ( specs.contains( "fieldReferences" ) )
    {
        if ( !specs.at( "fieldReferences" ).is_array() )
            throw std::invalid_argument( "fieldReferences must be an array" );
        for ( auto const& fieldSpec : specs.at( "fieldReferences" ) )
            config.fieldReferences.push_back( parseFieldReference( fieldSpec, "fieldReferences[]" ) );
    }

    return config;
}

template<int Dim>
node_type
toNode( point_type<Dim> const& point )
{
    node_type node( Dim );
    for ( int i = 0; i < Dim; ++i )
        node( i ) = point( i );
    return node;
}

template<int Dim, typename SpaceType, typename ElementType>
double
evaluateDisplacementAtPoint( SpaceType const& Vh,
                             ElementType const& u,
                             point_type<Dim> const& point,
                             point_type<Dim> const& direction )
{
    auto ctx = Vh->context();
    ctx.add( toNode( point ) );
    auto const value = u.evaluate( ctx, false );
    if ( value.size() < Dim )
        throw std::runtime_error( "displacement evaluation returned fewer components than expected" );

    double computed = 0.0;
    for ( int i = 0; i < Dim; ++i )
        computed += direction( i ) * value( i );
    return computed;
}

template<int Dim, typename SpaceType, typename ElementType, typename RangeType>
double
evaluateMarkedPointComponent( SpaceType const& Vh,
                              ElementType const& u,
                              RangeType const& tipRange,
                              std::string const& tipMarker,
                              int component )
{
    size_type nTipPoints = nelements( tipRange, true );
    if ( nTipPoints == 0 )
        throw std::invalid_argument( "cantilever tip marker '" + tipMarker + "' selects no point" );

    auto tipComponentDofs = Vh->dofs( tipRange, static_cast<ComponentType>( component ), false );
    double tipDisplacementSum = 0.0;
    size_type nTipComponentDofs = 0;
    for ( auto dofId : tipComponentDofs )
    {
        if ( Vh->dof()->dofGlobalProcessIsGhost( dofId ) )
            continue;
        tipDisplacementSum += u( dofId );
        ++nTipComponentDofs;
    }
    mpi::all_reduce( Environment::worldComm(), mpi::inplace( tipDisplacementSum ), std::plus<double>() );
    mpi::all_reduce( Environment::worldComm(), mpi::inplace( nTipComponentDofs ), std::plus<size_type>() );
    if ( nTipComponentDofs != nTipPoints )
        throw std::runtime_error( "cantilever tip marker '" + tipMarker + "' resolved to " +
                                  std::to_string( nTipComponentDofs ) + " active component dofs but " +
                                  std::to_string( nTipPoints ) + " marked points" );
    return tipDisplacementSum/nTipComponentDofs;
}

template<typename ValuesType>
voigt_type
toVoigtFromSymmetricStorage( ValuesType const& values, bool engineeringShear )
{
    double const shearScale = engineeringShear ? 2.0 : 1.0;
    voigt_type result = voigt_type::Zero();

    if ( values.size() >= 9 )
    {
        result << values( 0 ),                                      // xx
                  values( 4 ),                                      // yy
                  values( 8 ),                                      // zz
                  shearScale*0.5*( values( 1 ) + values( 3 ) ),     // xy
                  shearScale*0.5*( values( 2 ) + values( 6 ) ),     // xz
                  shearScale*0.5*( values( 5 ) + values( 7 ) );     // yz
        return result;
    }

    if ( values.size() < 6 )
        throw std::runtime_error( "expected at least six symmetric tensor components" );
    result << values( 0 ),                  // xx
              values( 3 ),                  // yy
              values( 5 ),                  // zz
              shearScale*values( 1 ),       // xy
              shearScale*values( 2 ),       // xz
              shearScale*values( 4 );       // yz
    return result;
}

inline bool
reportScalarComparison( std::string const& label,
                        double computed,
                        double reference,
                        Tolerance const& tolerance,
                        bool checkEnabled,
                        std::ostream& os )
{
    double const absoluteError = computed - reference;
    double const relError = relativeError( computed, reference );
    bool ok = true;
    if ( checkEnabled && tolerance.hasAbsolute )
        ok = ok && std::abs( absoluteError ) <= tolerance.absolute;
    if ( checkEnabled && tolerance.hasRelative )
        ok = ok && relError <= tolerance.relative;

    os << label << " computed = " << formatDouble( computed ) << "\n"
       << label << " reference = " << formatDouble( reference ) << "\n"
       << label << " absolute error = " << formatDouble( absoluteError ) << "\n"
       << label << " relative error = " << formatDouble( relError ) << "\n";
    if ( tolerance.hasAbsolute )
        os << label << " absolute tolerance = " << formatDouble( tolerance.absolute ) << "\n";
    if ( tolerance.hasRelative )
        os << label << " relative tolerance = " << formatDouble( tolerance.relative ) << "\n";
    if ( checkEnabled && ( tolerance.hasAbsolute || tolerance.hasRelative ) )
        os << label << " reference check = " << ( ok ? "PASS" : "FAIL" ) << "\n";
    else if ( !checkEnabled )
        os << label << " reference check = REPORT-ONLY (check disabled)\n";
    else
        os << label << " reference check = REPORT-ONLY (no tolerance configured)\n";
    return ok;
}

inline bool
reportScalarComparisonWithCombinedTolerance( std::string const& label,
                                             double computed,
                                             double reference,
                                             Tolerance const& tolerance,
                                             bool checkEnabled,
                                             std::ostream& os )
{
    double const absoluteError = computed - reference;
    double const relError = relativeError( computed, reference );
    bool const hasTolerance = tolerance.hasAbsolute || tolerance.hasRelative;
    double allowedError = 0.0;
    if ( tolerance.hasAbsolute )
        allowedError = std::max( allowedError, tolerance.absolute );
    if ( tolerance.hasRelative )
        allowedError = std::max( allowedError, tolerance.relative*std::max( std::abs( reference ), 1e-30 ) );
    bool const ok = !checkEnabled || !hasTolerance || std::abs( absoluteError ) <= allowedError;

    os << label << " computed = " << formatDouble( computed ) << "\n"
       << label << " reference = " << formatDouble( reference ) << "\n"
       << label << " absolute error = " << formatDouble( absoluteError ) << "\n"
       << label << " relative error = " << formatDouble( relError ) << "\n";
    if ( tolerance.hasAbsolute )
        os << label << " absolute tolerance = " << formatDouble( tolerance.absolute ) << "\n";
    if ( tolerance.hasRelative )
        os << label << " relative tolerance = " << formatDouble( tolerance.relative ) << "\n";
    if ( hasTolerance )
        os << label << " effective absolute tolerance = " << formatDouble( allowedError ) << "\n";
    if ( checkEnabled && hasTolerance )
        os << label << " reference check = " << ( ok ? "PASS" : "FAIL" ) << "\n";
    else if ( !checkEnabled )
        os << label << " reference check = REPORT-ONLY (check disabled)\n";
    else
        os << label << " reference check = REPORT-ONLY (no tolerance configured)\n";
    return ok;
}

inline bool
reportVoigtComparison( std::string const& label,
                       voigt_type const& computed,
                       voigt_type const& reference,
                       Tolerance const& tolerance,
                       bool checkEnabled,
                       std::ostream& os )
{
    static constexpr std::array<char const*, 6> components = { "xx", "yy", "zz", "xy", "xz", "yz" };

    double maxAbs = 0.0;
    double maxRel = 0.0;
    bool ok = true;
    os << label << " computed = {";
    for ( int i = 0; i < 6; ++i )
        os << ( i ? "," : "" ) << formatDouble( computed( i ) );
    os << "}\n" << label << " reference = {";
    for ( int i = 0; i < 6; ++i )
        os << ( i ? "," : "" ) << formatDouble( reference( i ) );
    os << "}\n";

    for ( int i = 0; i < 6; ++i )
    {
        double const absoluteError = computed( i ) - reference( i );
        double const relError = relativeError( computed( i ), reference( i ) );
        maxAbs = std::max( maxAbs, std::abs( absoluteError ) );
        maxRel = std::max( maxRel, relError );
        if ( checkEnabled && tolerance.hasAbsolute )
            ok = ok && std::abs( absoluteError ) <= tolerance.absolute;
        if ( checkEnabled && tolerance.hasRelative )
            ok = ok && relError <= tolerance.relative;

        os << "  " << label << "." << components[i]
           << ": computed=" << formatDouble( computed( i ) )
           << ", reference=" << formatDouble( reference( i ) )
           << ", abs-error=" << formatDouble( absoluteError )
           << ", rel-error=" << formatDouble( relError ) << "\n";
    }

    os << label << " max absolute error = " << formatDouble( maxAbs ) << "\n"
       << label << " max relative error = " << formatDouble( maxRel ) << "\n";
    if ( tolerance.hasAbsolute )
        os << label << " absolute tolerance = " << formatDouble( tolerance.absolute ) << "\n";
    if ( tolerance.hasRelative )
        os << label << " relative tolerance = " << formatDouble( tolerance.relative ) << "\n";
    if ( checkEnabled && ( tolerance.hasAbsolute || tolerance.hasRelative ) )
        os << label << " reference check = " << ( ok ? "PASS" : "FAIL" ) << "\n";
    else if ( !checkEnabled )
        os << label << " reference check = REPORT-ONLY (check disabled)\n";
    else
        os << label << " reference check = REPORT-ONLY (no tolerance configured)\n";
    return ok;
}

template<int Dim, typename SpaceType, typename ElementType>
int
checkProbeReference( Config<Dim> const& config,
                     SpaceType const& Vh,
                     ElementType const& u,
                     std::ostream& os = std::cout )
{
    if ( !config.hasProbeReference() )
        return 0;

    for ( auto const& probe : config.probes )
    {
        if ( !probe.hasProbe || !probe.hasReference )
            continue;

        double const computed = evaluateDisplacementAtPoint( Vh, u, probe.point, probe.direction );

        os << "probe reference check: " << probe.label << "\n";
        if ( !probe.description.empty() )
            os << "probe reference description = " << probe.description << "\n";
        bool const ok = reportScalarComparison( "probe displacement", computed, probe.value,
                                                probe.tolerance, probe.check, os );
        if ( !ok )
            throw std::runtime_error( "probe displacement reference check failed" );
    }
    return 0;
}

template<int Dim, typename SpaceType, typename ElementType>
Checker&
addProbeReferenceCheck( Checker& checker,
                        Config<Dim> const& config,
                        bool enabled,
                        SpaceType const& Vh,
                        ElementType const& u )
{
    return checker.add( "probe displacement reference", enabled && config.hasProbeReference(), [&config,&Vh,&u]() {
        return checkProbeReference( config, Vh, u );
    } );
}

inline double
cantileverEulerBernoulliReference( double youngModulus,
                                   CantileverReference<1> const& cantilever )
{
    if ( youngModulus <= 0.0 )
        throw std::invalid_argument( "cantilever Young modulus must be positive" );
    if ( cantilever.length <= 0.0 )
        throw std::invalid_argument( "cantilever.length must be positive" );
    if ( cantilever.height <= 0.0 )
        throw std::invalid_argument( "cantilever.height must be positive" );
    if ( cantilever.thickness <= 0.0 )
        throw std::invalid_argument( "cantilever.thickness must be positive" );

    double const inertia = cantilever.thickness*std::pow( cantilever.height, 3 )/12.0;
    return cantilever.tipForce*std::pow( cantilever.length, 3 )/( 3.0*youngModulus*inertia )
           + cantilever.tipMoment*std::pow( cantilever.length, 2 )/( 2.0*youngModulus*inertia );
}

template<int Dim>
double
cantileverEulerBernoulliReference( double youngModulus,
                                   CantileverReference<Dim> const& cantilever )
{
    CantileverReference<1> scalarCantilever;
    scalarCantilever.length = cantilever.length;
    scalarCantilever.height = cantilever.height;
    scalarCantilever.thickness = cantilever.thickness;
    scalarCantilever.tipForce = cantilever.tipForce;
    scalarCantilever.tipMoment = cantilever.tipMoment;
    return cantileverEulerBernoulliReference( youngModulus, scalarCantilever );
}

template<int Dim, typename SpaceType, typename ElementType, typename MeshType>
int
checkCantileverReferences( Config<Dim> const& config,
                           SpaceType const& Vh,
                           ElementType const& u,
                           MeshType const& mesh,
                           double youngModulus,
                           std::string const& target = "hexa8",
                           std::ostream& os = std::cout )
{
    for ( auto const& cantilever : config.cantilevers )
    {
        if ( cantilever.component < 0 || cantilever.component >= Dim )
            throw std::invalid_argument( "cantilever.component must be in [0," + std::to_string( Dim-1 ) + "]" );

        double computed = 0.0;
        if ( cantilever.hasPoint )
            computed = evaluateDisplacementAtPoint( Vh, u, cantilever.point, cantilever.direction );
        else
            computed = evaluateMarkedPointComponent<Dim>( Vh, u, markedpoints( mesh, cantilever.tipMarker ),
                                                          cantilever.tipMarker, cantilever.component );

        double const reference = cantileverEulerBernoulliReference( youngModulus, cantilever );
        bool checkEnabled = cantilever.check;
        Tolerance tolerance = cantilever.tolerance;
        if ( auto targetIt = cantilever.targets.find( target ); targetIt != cantilever.targets.end() )
        {
            auto const& options = targetIt->second;
            if ( options.hasCheck )
                checkEnabled = options.check;
            if ( hasTolerance( options.tolerance ) )
                tolerance = options.tolerance;
        }

        os << "Euler-Bernoulli cantilever reference check: " << cantilever.name << "\n";
        if ( !cantilever.description.empty() )
            os << "Euler-Bernoulli cantilever description = " << cantilever.description << "\n";
        os << "Euler-Bernoulli cantilever target = " << target << "\n";
        os << "Euler-Bernoulli cantilever length = " << formatDouble( cantilever.length ) << "\n"
           << "Euler-Bernoulli cantilever height = " << formatDouble( cantilever.height ) << "\n"
           << "Euler-Bernoulli cantilever thickness = " << formatDouble( cantilever.thickness ) << "\n"
           << "Euler-Bernoulli cantilever tip force = " << formatDouble( cantilever.tipForce ) << "\n"
           << "Euler-Bernoulli cantilever tip moment = " << formatDouble( cantilever.tipMoment ) << "\n";

        bool const ok = reportScalarComparisonWithCombinedTolerance( "cantilever tip displacement", computed, reference,
                                                                     tolerance, checkEnabled, os );
        if ( !ok )
            throw std::runtime_error( "Euler-Bernoulli cantilever reference check failed" );
    }
    return 0;
}

template<int Dim, typename SpaceType, typename ElementType, typename MeshType>
Checker&
addCantileverReferenceChecks( Checker& checker,
                              Config<Dim> const& config,
                              bool enabled,
                              SpaceType const& Vh,
                              ElementType const& u,
                              MeshType const& mesh,
                              double youngModulus,
                              std::string target = "hexa8" )
{
    return checker.add( "Euler-Bernoulli cantilever", enabled && config.hasCantileverReferences(),
                        [&config,&Vh,&u,&mesh,youngModulus,target=std::move( target )]() {
        return checkCantileverReferences( config, Vh, u, mesh, youngModulus, target );
    } );
}

template<int Dim, typename SpaceType, typename FieldElementType>
voigt_type
evaluateSymmetricFieldReference( SpaceType const& Sh,
                                 FieldElementType const& field,
                                 point_type<Dim> const& point,
                                 bool engineeringShear )
{
    auto componentValue = [&field,&point]( ComponentType c1, ComponentType c2 )
    {
        auto component = field.comp( c1, c2 );
        auto ctx = component.functionSpace()->context();
        ctx.add( toNode( point ) );
        auto values = component.evaluate( ctx, true );
        if ( values.size() < 1 )
            throw std::runtime_error( "expected at least one scalar component value" );
        return values( 0 );
    };

    double const shearScale = engineeringShear ? 2.0 : 1.0;
    voigt_type result = voigt_type::Zero();
    result << componentValue( ComponentType::X, ComponentType::X ),
              componentValue( ComponentType::Y, ComponentType::Y ),
              componentValue( ComponentType::Z, ComponentType::Z ),
              shearScale*componentValue( ComponentType::X, ComponentType::Y ),
              shearScale*componentValue( ComponentType::X, ComponentType::Z ),
              shearScale*componentValue( ComponentType::Y, ComponentType::Z );
    return result;
}

template<int Dim, typename SpaceType, typename EpsilonElementType, typename SigmaElementType>
int
checkFieldReferences( Config<Dim> const& config,
                      SpaceType const& Sh,
                      EpsilonElementType const& epsilonh,
                      SigmaElementType const& sigmah,
                      std::string const& target = "hexa8",
                      std::ostream& os = std::cout )
{
    for ( auto const& fieldReference : config.fieldReferences )
    {
        auto const targetIt = fieldReference.targets.find( target );
        if ( targetIt == fieldReference.targets.end() )
            throw std::invalid_argument( "field reference '" + fieldReference.name + "' has no target '" + target + "'" );

        auto const computedEpsilon = evaluateSymmetricFieldReference( Sh, epsilonh, fieldReference.point, true );
        auto const computedSigma = evaluateSymmetricFieldReference( Sh, sigmah, fieldReference.point, false );

        os << "field reference check '" << fieldReference.name << "' target '" << target << "'\n";
        if ( !fieldReference.description.empty() )
            os << "field reference description = " << fieldReference.description << "\n";

        bool ok = true;
        auto const& reference = targetIt->second;
        if ( reference.hasEpsilon )
            ok = reportVoigtComparison( "epsilon", computedEpsilon, reference.epsilon,
                                        fieldReference.tolerance, fieldReference.check, os ) && ok;
        if ( reference.hasSigma )
            ok = reportVoigtComparison( "sigma", computedSigma, reference.sigma,
                                        fieldReference.tolerance, fieldReference.check, os ) && ok;
        if ( !ok )
            throw std::runtime_error( "field reference '" + fieldReference.name + "' check failed" );
    }
    return 0;
}

template<int Dim, typename SpaceType, typename EpsilonElementType, typename SigmaElementType>
Checker&
addFieldReferenceChecks( Checker& checker,
                         Config<Dim> const& config,
                         bool enabled,
                         SpaceType const& Sh,
                         EpsilonElementType const& epsilonh,
                         SigmaElementType const& sigmah,
                         std::string target = "hexa8" )
{
    return checker.add( "field reference values", enabled && config.hasFieldReferences(),
                        [&config,&Sh,&epsilonh,&sigmah,target=std::move( target )]() {
        return checkFieldReferences( config, Sh, epsilonh, sigmah, target );
    } );
}

template<int Dim>
class ElasticityReferenceChecker : public Checker
{
public:
    using Checker::Checker;

    ElasticityReferenceChecker& load( nl::json const& specs )
    {
        return this->setConfig( jsonConfig<Dim>( specs ) );
    }

    ElasticityReferenceChecker& setConfig( Config<Dim> config )
    {
        M_config = std::move( config );
        M_target = M_config.target;
        return *this;
    }

    Config<Dim> const& config() const { return M_config; }

    ElasticityReferenceChecker& setTarget( std::string target )
    {
        M_target = std::move( target );
        return *this;
    }

    std::string const& target() const { return M_target; }

    template<typename SpaceType, typename ElementType>
    ElasticityReferenceChecker& addDisplacementProbe( SpaceType const& Vh,
                                                      ElementType const& u,
                                                      bool enabled = true )
    {
        this->add( "probe displacement reference", enabled && M_config.hasProbeReference(), [this,&Vh,&u]() {
            return checkProbeReference( M_config, Vh, u );
        } );
        return *this;
    }

    template<typename MeshType, typename SpaceType, typename ElementType>
    ElasticityReferenceChecker& addCantileverChecks( MeshType const& mesh,
                                                     SpaceType const& Vh,
                                                     ElementType const& u,
                                                     double youngModulus,
                                                     bool enabled = true )
    {
        this->add( "Euler-Bernoulli cantilever", enabled && M_config.hasCantileverReferences(),
                   [this,&mesh,&Vh,&u,youngModulus]() {
            return checkCantileverReferences( M_config, Vh, u, mesh, youngModulus, M_target );
        } );
        return *this;
    }

    template<typename MeshType, typename ElementType>
    ElasticityReferenceChecker& addSmallStrainStressChecks( MeshType const& mesh,
                                                            ElementType const& u,
                                                            double lambda,
                                                            double mu,
                                                            bool enabled = true )
    {
        if constexpr ( Dim == 3 )
        {
            this->add( "field reference values", enabled && M_config.hasFieldReferences(), [this,&mesh,&u,lambda,mu]() {
                auto Sh = Pdhms<0>( mesh );
                auto epsilonTensor = cst( 0.5 )*( gradv( u ) + trans( gradv( u ) ) );
                auto sigmaTensor = cst( lambda )*trace( epsilonTensor )*eye<3,3>() + cst( 2.0*mu )*epsilonTensor;
                auto epsilonh = vf::project( _space=Sh, _range=elements( mesh ), _expr=epsilonTensor );
                auto sigmah = vf::project( _space=Sh, _range=elements( mesh ), _expr=sigmaTensor );
                return checkFieldReferences( M_config, Sh, epsilonh, sigmah, M_target );
            } );
        }
        else
        {
            this->add( "field reference values", enabled && M_config.hasFieldReferences(), []() -> int {
                throw std::invalid_argument( "fieldReferences are only supported by qs_elasticity in 3D" );
            } );
        }
        return *this;
    }

private:
    Config<Dim> M_config;
    std::string M_target = "hexa8";
};

} // namespace Feel::Quickstart::ElasticityChecks

#endif
